// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataHandling/LoadBankFromDiskTask.h"
#include "MantidDataHandling/BankPulseTimes.h"
#include "MantidDataHandling/DefaultEventLoader.h"
#include "MantidDataHandling/LoadEventNexus.h"
#include "MantidDataHandling/ProcessBankCompressed.h"
#include "MantidDataHandling/ProcessBankData.h"
#include "MantidDataHandling/PulseIndexer.h"
#include "MantidKernel/ParallelMinMax.h"
#include "MantidKernel/Timer.h"
#include "MantidKernel/Unit.h"
#include "MantidKernel/VectorHelper.h"
#include "MantidNexus/NexusException.h"
#include "MantidNexus/NexusFile.h"
#include "MantidNexus/NexusIOHelper.h"

#include "tbb/blocked_range.h"
#include "tbb/global_control.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>

namespace {
// this is used for unit conversion to correct units
const std::string MICROSEC("microseconds");

/** EventsPerDetIdCounter
 * A functor for use with tbb::parallel_reduce, which counts the events on each detector ID in [minId, maxId] over a
 * subrange of the event ids. Each split counts into its own vector, and these are then added together.
 */
class EventsPerDetIdCounter {
  std::vector<uint32_t> const *m_detIds;
  uint32_t m_minId;
  uint32_t m_maxId;

public:
  std::vector<size_t> counts;

  EventsPerDetIdCounter(std::vector<uint32_t> const *detIds, const uint32_t minId, const uint32_t maxId)
      : m_detIds(detIds), m_minId(minId), m_maxId(maxId), counts(static_cast<size_t>(maxId - minId) + 1, 0) {}

  // start the split with its own empty counts
  EventsPerDetIdCounter(EventsPerDetIdCounter &other, tbb::split)
      : m_detIds(other.m_detIds), m_minId(other.m_minId), m_maxId(other.m_maxId), counts(other.counts.size(), 0) {}

  void operator()(tbb::blocked_range<size_t> const &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      const auto detId = (*m_detIds)[i];
      if (detId >= m_minId && detId <= m_maxId)
        ++counts[detId - m_minId];
    }
  }

  void join(EventsPerDetIdCounter const &other) {
    std::transform(counts.cbegin(), counts.cend(), other.counts.cbegin(), counts.begin(), std::plus<size_t>());
  }
};

/** Count the events on each detector ID in [minId, maxId], using the first numEvents entries of detIds
 * @return the counts, where index 0 is minId
 */
std::vector<size_t> countEventsPerDetId(std::vector<uint32_t> const &detIds, const size_t numEvents,
                                        const uint32_t minId, const uint32_t maxId) {
  // large grains keep the number of splits, each with its own counts, small
  constexpr size_t GRAINSIZE{1 << 20};
  EventsPerDetIdCounter counter(&detIds, minId, maxId);
  tbb::parallel_reduce(tbb::blocked_range<size_t>(0, std::min(numEvents, detIds.size()), GRAINSIZE), counter);
  return std::move(counter.counts);
}

/** Find the detector ID to split [minId, maxId] at so each side has about half the events
 * @param counts :: events per detector ID, where index 0 is minId
 * @param minId :: first detector ID in counts
 * @return the last detector ID of the first half
 */
uint32_t balancedMidId(std::vector<size_t> const &counts, const uint32_t minId) {
  const size_t total = std::accumulate(counts.cbegin(), counts.cend(), static_cast<size_t>(0));
  const size_t half = (total + 1) / 2;
  size_t cumulative = 0;
  for (size_t i = 0; i < counts.size(); ++i) {
    cumulative += counts[i];
    if (cumulative >= half)
      return minId + static_cast<uint32_t>(i);
  }
  return minId + static_cast<uint32_t>(counts.size() - 1);
}

/// Events in each range that fillInPlace fills in parallel, so banks with fewer events are filled in one range
constexpr size_t EVENTS_PER_FILL_RANGE{size_t(1) << 24};

/** How many contiguous ranges of events to fill a bank in, in parallel: one per EVENTS_PER_FILL_RANGE events, but no
 * more than the cores TBB may use (which follows MultiThreaded.MaxCores), and no more than numEvents / numDetIds. Each
 * range keeps a count and a write position for every detector ID, so the last limit keeps those smaller than the
 * events the range writes.
 * @param numEvents :: events in the bank
 * @param numDetIds :: detector IDs the bank's events are filled into
 */
size_t numFillRanges(const size_t numEvents, const size_t numDetIds) {
  const size_t byEvents = (numEvents + EVENTS_PER_FILL_RANGE - 1) / EVENTS_PER_FILL_RANGE;
  const size_t byCores = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
  const size_t byMemory = numEvents / std::max<size_t>(numDetIds, 1);
  return std::max<size_t>(1, std::min({byEvents, byCores, byMemory}));
}

/** The pulses with events in a bank's arrays, split into contiguous ranges holding about the same number of events.
 *
 * A PulseIndexer over part of the arrays does not find the right pulses, so the split is made between the pulses of
 * one PulseIndexer over all of them.
 */
class PulseRanges {
public:
  PulseRanges(const Mantid::DataHandling::PulseIndexer &indexer, const size_t numRanges)
      : m_rangeStart(numRanges + 1, 0) {
    size_t eventsInPulses = 0;
    for (const auto &pulse : indexer) {
      if (pulse.eventIndexStop > pulse.eventIndexStart) {
        m_pulses.push_back(pulse);
        eventsInPulses += pulse.eventIndexStop - pulse.eventIndexStart;
      }
    }
    std::fill(m_rangeStart.begin() + 1, m_rangeStart.end(), m_pulses.size());
    size_t eventsSoFar = 0;
    size_t range = 1;
    for (size_t k = 0; k < m_pulses.size() && range < numRanges; ++k) {
      eventsSoFar += m_pulses[k].eventIndexStop - m_pulses[k].eventIndexStart;
      while (range < numRanges && eventsSoFar * numRanges >= eventsInPulses * range)
        m_rangeStart[range++] = k + 1;
    }
  }

  size_t size() const { return m_rangeStart.size() - 1; }

  /// Call func(arrayIndex, pulseTime) for each event in one range, in order
  template <typename Func>
  void forEachEvent(const size_t range, const Mantid::DataHandling::BankPulseTimes &pulseTimes, Func &&func) const {
    for (size_t k = m_rangeStart[range]; k < m_rangeStart[range + 1]; ++k) {
      const auto &pulse = m_pulses[k];
      const auto &pulseTime = pulseTimes.pulseTime(pulse.pulseIndex);
      for (size_t i = pulse.eventIndexStart; i < pulse.eventIndexStop; ++i)
        func(i, pulseTime);
    }
  }

private:
  std::vector<Mantid::DataHandling::PulseIndexer::IteratorValue> m_pulses;
  /// range r has pulses [m_rangeStart[r], m_rangeStart[r + 1])
  std::vector<size_t> m_rangeStart;
};

using EventLists = std::vector<std::vector<Mantid::Types::Event::TofEvent> *>;
using Mantid::Types::Event::TofEvent;

/** Count each range's events for each list, in parallel
 * @param isKept :: isKept(id, tof) is whether the filters keep an event with an ID in [minId, minId + width)
 * @return counts[range][id - minId]
 */
template <typename IsKept>
std::vector<std::vector<size_t>>
countEventsPerRange(const PulseRanges &ranges, const Mantid::DataHandling::BankPulseTimes &pulseTimes,
                    std::vector<uint32_t> const &ids, std::vector<float> const &tofs, const IsKept &isKept,
                    const EventLists &lists, const uint32_t minId, const size_t width) {
  std::vector<std::vector<size_t>> rangeCounts(ranges.size(), std::vector<size_t>(width, 0));
  tbb::parallel_for(size_t(0), ranges.size(), [&](const size_t range) {
    auto &counts = rangeCounts[range];
    ranges.forEachEvent(range, pulseTimes, [&](const size_t i, const auto &) {
      const uint32_t id = ids[i];
      if (isKept(id, static_cast<double>(tofs[i])) && lists[id])
        ++counts[id - minId];
    });
  });
  return rangeCounts;
}

/// Whether two detector IDs with events, where index 0 of totals is minId, share one list
bool listsAreShared(const EventLists &lists, const uint32_t minId, std::vector<size_t> const &totals) {
  std::vector<const void *> used;
  for (size_t i = 0; i < totals.size(); ++i) {
    if (totals[i] > 0)
      used.emplace_back(lists[minId + i]);
  }
  std::sort(used.begin(), used.end());
  return std::adjacent_find(used.cbegin(), used.cend()) != used.cend();
}

/** Grow each list once to fit its new events, in parallel, and give every range its start in each list, after any
 * events already there and the earlier ranges' events
 * @return cursors[range][id - minId], null for detector IDs with no events
 */
std::vector<std::vector<TofEvent *>> growLists(const EventLists &lists, const uint32_t minId,
                                               std::vector<size_t> const &totals,
                                               std::vector<std::vector<size_t>> const &rangeCounts) {
  std::vector<std::vector<TofEvent *>> rangeCursors(rangeCounts.size(),
                                                    std::vector<TofEvent *>(totals.size(), nullptr));
  tbb::parallel_for(size_t(0), totals.size(), [&](const size_t i) {
    if (totals[i] == 0)
      return;
    auto *list = lists[minId + i];
    const size_t existing = list->size();
    list->resize(existing + totals[i]);
    auto *position = list->data() + existing;
    for (size_t range = 0; range < rangeCounts.size(); ++range) {
      rangeCursors[range][i] = position;
      position += rangeCounts[range][i];
    }
  });
  return rangeCursors;
}

/** Write each range's events from its own start in every list, in parallel
 * @return the time-of-flight limits and counts of all the ranges
 */
template <typename IsKept>
Mantid::DataHandling::TofStats
fillRanges(const PulseRanges &ranges, const Mantid::DataHandling::BankPulseTimes &pulseTimes,
           std::vector<uint32_t> const &ids, std::vector<float> const &tofs, const IsKept &isKept, const uint32_t minId,
           std::vector<std::vector<TofEvent *>> &rangeCursors) {
  std::vector<Mantid::DataHandling::TofStats> stats(ranges.size());
  tbb::parallel_for(size_t(0), ranges.size(), [&](const size_t range) {
    auto &cursors = rangeCursors[range];
    auto &stat = stats[range];
    ranges.forEachEvent(range, pulseTimes, [&](const size_t i, const auto &pulseTime) {
      const uint32_t id = ids[i];
      const auto tof = static_cast<double>(tofs[i]);
      if (!isKept(id, tof))
        return;
      auto *&cursor = cursors[id - minId];
      if (cursor) {
        *cursor = TofEvent(tof, pulseTime);
        ++cursor;
      } else {
        ++stat.discardedEvents;
      }
      stat.add(tof);
    });
  });
  Mantid::DataHandling::TofStats total;
  for (const auto &stat : stats)
    total.merge(stat);
  return total;
}

} // namespace

namespace Mantid::DataHandling {

/** Constructor
 *
 * @param loader :: Handle to the main loader
 * @param entry_name :: The address of the bank to load
 * @param entry_type :: The classtype of the entry to load
 * @param numEvents :: The number of events in the bank.
 * @param oldNeXusFileNames :: Identify if file is of old variety.
 * @param prog :: an optional Progress object
 * @param scheduler :: the ThreadScheduler that runs this task.
 * @param framePeriodNumbers :: Period numbers corresponding to each frame
 */
LoadBankFromDiskTask::LoadBankFromDiskTask(DefaultEventLoader &loader, std::string entry_name, std::string entry_type,
                                           const std::size_t numEvents, const bool oldNeXusFileNames,
                                           API::Progress *prog, Kernel::ThreadScheduler &scheduler,
                                           std::vector<int> framePeriodNumbers)
    : m_loader(loader), entry_name(std::move(entry_name)), entry_type(std::move(entry_type)), prog(prog),
      scheduler(scheduler), m_loadError(false), m_have_weight(false),
      m_framePeriodNumbers(std::move(framePeriodNumbers)) {
  m_cost = static_cast<double>(numEvents);

  // some field names changed over time
  m_timeOfFlightFieldName = oldNeXusFileNames ? "event_time_of_flight" : "event_time_offset";
  m_detIdFieldName = oldNeXusFileNames ? "event_pixel_id" : "event_id";

  // detector id range
  m_min_id = std::numeric_limits<uint32_t>::max();
  m_max_id = 0;
}

/** Load the pulse times, if needed. This sets thisBankPulseTimes to the right pointer.
 */
void LoadBankFromDiskTask::loadPulseTimes(Nexus::File &file) {
  try {
    // First, get info about the event_time_zero field in this bank
    file.openData("event_time_zero");
  } catch (const Nexus::Exception &) {
    // Field not found error is most likely.
    // Use the "proton_charge" das logs.
    thisBankPulseTimes = m_loader.alg->m_allBanksPulseTimes;
    return;
  }
  std::string thisStartTime;
  size_t thispulseTimes = 0;
  // If the offset is not present, use Unix epoch
  if (!file.hasAttr("offset")) {
    thisStartTime = "1970-01-01T00:00:00Z";
    m_loader.alg->getLogger().warning() << "In loadPulseTimes: no ISO8601 offset attribute provided for "
                                           "event_time_zero, using UNIX epoch instead\n";
  } else {
    file.getAttr("offset", thisStartTime);
  }

  if (!file.getInfo().dims.empty())
    thispulseTimes = static_cast<size_t>(file.getInfo().dims[0]);
  file.closeData();

  // Now, we look through existing ones to see if it is already loaded. Other bank tasks search and add to the list
  // at the same time, so it is only touched under its mutex.
  const auto findLoaded = [this, thispulseTimes, &thisStartTime]() {
    const auto &loaded = m_loader.m_bankPulseTimes;
    const auto found = std::find_if(loaded.cbegin(), loaded.cend(), [&](const auto &bankPulseTime) {
      return bankPulseTime->equals(thispulseTimes, thisStartTime);
    });
    return found == loaded.cend() ? nullptr : *found;
  };
  {
    std::lock_guard<std::mutex> lock(m_loader.m_bankPulseTimesMutex);
    if (auto loaded = findLoaded()) {
      thisBankPulseTimes = std::move(loaded);
      return;
    }
  }

  // Not found? Need to load and add it. Load outside the lock so other banks can read at the same time, then use
  // whichever copy is in the list, in case another bank added the same pulse times meanwhile.
  auto pulseTimes = std::make_shared<BankPulseTimes>(file, m_framePeriodNumbers);
  std::lock_guard<std::mutex> lock(m_loader.m_bankPulseTimesMutex);
  if (auto loaded = findLoaded()) {
    thisBankPulseTimes = std::move(loaded);
  } else {
    m_loader.m_bankPulseTimes.emplace_back(pulseTimes);
    thisBankPulseTimes = std::move(pulseTimes);
  }
}

/** Load the event_index field
 * (a list of size of # of pulses giving the index in the event list for that pulse)
 * @param file :: File handle for the NeXus file
 */
std::unique_ptr<std::vector<uint64_t>> LoadBankFromDiskTask::loadEventIndex(Nexus::File &file) {
  // Get the event_index (a list of size of # of pulses giving the index in
  // the event list for that pulse) as a uint64 vector.
  // The Nexus standard does not specify if this is to be 32-bit or 64-bit
  // integers, so we use the NeXusIOHelper to do the conversion on the fly.
  auto event_index =
      std::make_unique<std::vector<uint64_t>>(Nexus::IOHelper::readNexusVector<uint64_t>(file, "event_index"));

  // Look for the sign that the bank is empty
  if (event_index->size() == 1) {
    if (event_index->at(0) == 0) {
      // One entry, only zero. This means NO events in this bank.
      m_loadError = true;
      m_loader.alg->getLogger().debug() << "Bank " << entry_name << " is empty.\n";
    }
  }

  return event_index;
}

/** Open the event_id field and validate the contents
 *
 * @param file :: File handle for the NeXus file
 * @param start_event :: set to the index of the first event
 * @param stop_event :: set to the index of the last event + 1
 * @param start_event_index ::  (a list of size of # of pulses giving the index in
 * the event list for that pulse)
 */
void LoadBankFromDiskTask::prepareEventId(Nexus::File &file, uint64_t &start_event, uint64_t &stop_event,
                                          const uint64_t &start_event_index) {
  // Get the list of pixel ID's
  file.openData(m_detIdFieldName);

  // By default, use all available indices
  start_event = start_event_index;
  Nexus::Info id_info = file.getInfo();
  // dims[0] can be negative in ISIS meaning 2^32 + dims[0]. Take that into
  // account
  uint64_t dim0 = recalculateDataSize(id_info.dims[0]);
  stop_event = dim0;

  // We are loading part - work out the event number range
  if (m_loader.chunk != EMPTY_INT()) {
    start_event = (m_loader.chunk - m_loader.firstChunkForBank) * (m_loader.eventsPerChunk);
    // Don't change stop_event for the final chunk
    if (start_event + m_loader.eventsPerChunk < stop_event)
      stop_event = start_event + m_loader.eventsPerChunk;
  }

  // Make sure it is within range
  if (stop_event > dim0)
    stop_event = dim0;

  m_loader.alg->getLogger().debug() << entry_name << ": start_event " << start_event << " stop_event " << stop_event
                                    << "\n";
}

/** Load the event_id field, which has been opened
 * @param file An Nexus::File object opened at the correct group
 * @returns A new array containing the event Ids for this bank
 */
std::unique_ptr<std::vector<uint32_t>> LoadBankFromDiskTask::loadEventId(Nexus::File &file) {
  // This is the data size
  Nexus::Info id_info = file.getInfo();
  const Nexus::dimsize_t dim0 = recalculateDataSize(id_info.dims[0]);

  // Check that the required space is there in the file.
  if (dim0 < m_loadSize[0] + m_loadStart[0]) {
    m_loader.alg->getLogger().warning() << "Entry " << entry_name << "'s event_id field is too small (" << dim0
                                        << ") to load the desired data size (" << m_loadSize[0] + m_loadStart[0]
                                        << ").\n";
    m_loadError = true;
  }

  // Now we allocate the required arrays, sized for only the events being loaded
  auto event_id = std::make_unique<std::vector<uint32_t>>(m_loadError ? 0 : m_loadSize[0]);

  if (!m_loadError) {
    Nexus::IOHelper::readNexusSlab<uint32_t, Nexus::IOHelper::Narrowing::Prevent>(*event_id, file, m_detIdFieldName,
                                                                                  m_loadStart, m_loadSize);
    file.closeData();

    // determine the range of pixel ids
    {
      const auto [min_id, max_id] = Mantid::Kernel::parallel_minmax<uint32_t>(event_id);
      m_min_id = min_id;
      m_max_id = max_id;
    }

    if (m_min_id > static_cast<uint32_t>(m_loader.eventid_max)) {
      // All the detector IDs in the bank are higher than the highest 'known'
      // (from the IDF)
      // ID. Setting this will abort the loading of the bank.
      m_loadError = true;
    }
    // fixup the minimum pixel id in the case that it's lower than the lowest
    // 'known' id. We test this by checking that when we add the offset we
    // would not get a negative index into the vector. Note that m_min_id is
    // a uint so we have to be cautious about adding it to an int which may be
    // negative.
    if (static_cast<int32_t>(m_min_id) + m_loader.pixelID_to_wi_offset < 0) {
      m_min_id = static_cast<uint32_t>(abs(m_loader.pixelID_to_wi_offset));
    }
    // fixup the maximum pixel id in the case that it's higher than the
    // highest 'known' id
    if (m_max_id > static_cast<uint32_t>(m_loader.eventid_max))
      m_max_id = static_cast<uint32_t>(m_loader.eventid_max);
  }
  return event_id;
}

/** Open and load the times-of-flight data
 * @param file An Nexus::File object opened at the correct group
 * @returns A new array containing the time of flights for this bank
 */
std::unique_ptr<std::vector<float>> LoadBankFromDiskTask::loadTof(Nexus::File &file) {
  // Get the list of event_time_of_flight's
  file.openData(m_timeOfFlightFieldName);

  // This is the data size
  // Check that the required space is there in the file.
  Nexus::Info tof_info = file.getInfo();
  uint64_t tof_dim0 = recalculateDataSize(tof_info.dims[0]);
  if (tof_dim0 < m_loadSize[0] + m_loadStart[0]) {
    m_loader.alg->getLogger().warning() << "Entry " << entry_name
                                        << "'s event_time_offset field is too small "
                                           "to load the desired data.\n";
    m_loadError = true;
  }

  // Allocate the array, sized for only the events being loaded
  auto event_time_of_flight = std::make_unique<std::vector<float>>(m_loadError ? 0 : m_loadSize[0]);

  // Mantid assumes event_time_offset to be float.
  // Nexus only requires event_time_offset to be a NXNumber.
  // We thus have to consider 32-bit or 64-bit options, and we
  // explicitly allow downcasting using the additional AllowDowncasting
  // template argument.
  // the memory is allocated earlier in the function
  Nexus::IOHelper::readNexusSlab<float, Nexus::IOHelper::Narrowing::Allow>(
      *event_time_of_flight, file, m_timeOfFlightFieldName, m_loadStart, m_loadSize);
  std::string tof_unit;
  try {
    file.getAttr("units", tof_unit);
  } catch (Nexus::Exception const &) {
  }
  file.closeData();

  // Convert Tof to microseconds
  if (tof_unit != MICROSEC)
    Kernel::Units::timeConversionVector(*event_time_of_flight, tof_unit, MICROSEC);

  return event_time_of_flight;
}

/** Load weight of weigthed events if they exist
 * @param file An Nexus::File object opened at the correct group
 * @returns A new array containing the weights or a nullptr if the weights
 * are not present
 */
std::unique_ptr<std::vector<float>> LoadBankFromDiskTask::loadEventWeights(Nexus::File &file) {
  try {
    // First, get info about the event_weight field in this bank
    file.openData("event_weight");
  } catch (Nexus::Exception const &) {
    // Field not found error is most likely.
    m_have_weight = false;
    return std::unique_ptr<std::vector<float>>();
  }
  // OK, we've got them
  m_have_weight = true;

  // Allocate the array
  auto event_weight = std::make_unique<std::vector<float>>(m_loadSize[0]);

  Nexus::Info weight_info = file.getInfo();
  uint64_t weight_dim0 = recalculateDataSize(weight_info.dims[0]);
  if (weight_dim0 < m_loadSize[0] + m_loadStart[0]) {
    m_loader.alg->getLogger().warning() << "Entry " << entry_name
                                        << "'s event_weight field is too small to load the desired data.\n";
    m_loadError = true;
  }

  // Check that the type is what it is supposed to be
  if (weight_info.type == NXnumtype::FLOAT32)
    file.getSlab(event_weight->data(), m_loadStart, m_loadSize);
  else {
    m_loader.alg->getLogger().warning() << "Entry " << entry_name
                                        << "'s event_weight field is not FLOAT32! It will be skipped.\n";
    m_loadError = true;
  }

  if (!m_loadError) {
    file.closeData();
  }
  return event_weight;
}

void LoadBankFromDiskTask::run() {
  // timer for performance
  Mantid::Kernel::Timer timer;

  // These give the limits in each file as to which events we actually load
  // (when filtering by time).
  m_loadStart.resize(1, 0);
  m_loadSize.resize(1, 0);

  m_loadError = false;
  m_have_weight = m_loader.m_haveWeights;

  prog->report(entry_name + ": load from disk");

  // arrays to load into
  std::shared_ptr<std::vector<uint32_t>> event_id;
  std::shared_ptr<std::vector<float>> event_time_of_flight;
  std::shared_ptr<std::vector<float>> event_weight;
  std::shared_ptr<std::vector<uint64_t>> event_index;

  // Copy the algorithm's open file: shares the HDF5 handle so the descriptor is not rebuilt
  Nexus::File file(*m_loader.alg->m_file);
  try {
    // Navigate into the file
    file.openGroup(m_loader.alg->m_top_entry_name, "NXentry");
    // Open the bankN_event group
    file.openGroup(entry_name, entry_type);

    const bool needPulseInfo = (!m_loader.alg->compressEvents) || m_loader.alg->compressTolerance == 0 ||
                               m_loader.m_ws.nPeriods() > 1 || m_loader.alg->m_is_time_filtered ||
                               m_loader.alg->filter_bad_pulses || m_have_weight;

    // Load the event_index field.
    if (needPulseInfo)
      event_index = this->loadEventIndex(file);
    else
      event_index = nullptr;

    if (!m_loadError) {
      // Load and validate the pulse times
      if (needPulseInfo)
        this->loadPulseTimes(file);
      else
        thisBankPulseTimes = nullptr;
      // The event_index should be the same length as the pulse times from DAS
      // logs.
      if (event_index && event_index->size() != thisBankPulseTimes->numberOfPulses())
        m_loader.alg->getLogger().warning() << "Bank " << entry_name
                                            << " has a mismatch between the number of event_index entries "
                                               "and the number of pulse times in event_time_zero.\n";
      // Open and validate event_id field.
      uint64_t start_event = 0;
      uint64_t stop_event = 0;
      if (event_index)
        this->prepareEventId(file, start_event, stop_event, event_index->operator[](0));
      else
        this->prepareEventId(file, start_event, stop_event, 0);

      // These are the arguments to getSlab()
      m_loadStart[0] = start_event;
      m_loadSize[0] = stop_event - start_event;

      if ((m_loader.alg->compressEvents) || ((m_loadSize[0] > 0))) {
        if (m_loader.alg->getCancel()) {
          m_loader.alg->getLogger().error() << "Loading bank " << entry_name << " is cancelled.\n";
          m_loadError = true; // To allow cancelling the algorithm
        }

        // Load pixel IDs
        if (!m_loadError)
          event_id = this->loadEventId(file);

        // for compression the number of events needs to come from elsewhere
        if (!event_index)
          m_loadSize[0] = event_id->size();

        if (m_loader.alg->getCancel()) {
          m_loader.alg->getLogger().error() << "Loading bank " << entry_name << " is cancelled.\n";
          m_loadError = true; // To allow cancelling the algorithm
        }

        // And TOF.
        if (!m_loadError) {
          event_time_of_flight = this->loadTof(file);
          if (m_have_weight) {
            event_weight = this->loadEventWeights(file);
          }
        }
      } // Size is at least 1
      else {
        // Found a size that was 0 or less; stop processing
        m_loader.alg->getLogger().error()
            << "Loading bank " << entry_name << " is stopped due to either zero/negative loading size ("
            << m_loadStart[0] << ") or negative load start index (" << m_loadStart[0] << ")\n";
        m_loadError = true;
      }

    } // no error
  } // try block
  catch (std::exception &e) {
    m_loader.alg->getLogger().error() << "Error while loading bank " << entry_name << ":\n";
    m_loader.alg->getLogger().error() << e.what() << '\n';
    m_loadError = true;
  } catch (...) {
    m_loader.alg->getLogger().error() << "Unspecified error while loading bank " << entry_name << '\n';
    m_loadError = true;
  }

  // Close up the file even if errors occured.
  file.closeGroup();
  file.close();

  // Abort if anything failed
  if (m_loadError) {
    return;
  }

  const auto bank_size = m_max_id - m_min_id;
  const auto minSpectraToLoad = static_cast<uint32_t>(m_loader.alg->m_specMin);
  const auto maxSpectraToLoad = static_cast<uint32_t>(m_loader.alg->m_specMax);
  const auto emptyInt = static_cast<uint32_t>(EMPTY_INT());
  // check that if a range of spectra were requested that these fit within
  // this bank
  if (minSpectraToLoad != emptyInt && m_min_id < minSpectraToLoad) {
    if (minSpectraToLoad > m_max_id) { // the minimum spectra to load is more
                                       // than the max of this bank
      return;
    }
    // the min spectra to load is higher than the min for this bank
    m_min_id = minSpectraToLoad;
  }
  if (maxSpectraToLoad != emptyInt && m_max_id > maxSpectraToLoad) {
    if (maxSpectraToLoad < m_min_id) {
      // the maximum spectra to load is less than the minimum of this bank
      return;
    }
    // the max spectra to load is lower than the max for this bank
    m_max_id = maxSpectraToLoad;
  }
  if (m_min_id > m_max_id) {
    // the min is now larger than the max, this means the entire block of
    // spectra to load is outside this bank
    return;
  }

  // schedule the job to generate the event lists
  const bool useCompressed =
      (m_loader.alg->compressEvents) && (!event_weight) && (m_loader.alg->compressTolerance != 0);

  // Unweighted events from one period that are not compressed are written straight into lists sized to fit, in this
  // task. This needs the pre-count, and falls back to the processing tasks when it does not apply.
  const bool canFillInPlace = m_loader.precount && !m_have_weight && !m_loader.alg->compressEvents &&
                              m_loader.m_ws.nPeriods() == 1 && event_index;
  if (canFillInPlace &&
      fillInPlace(*event_id, *event_time_of_flight, static_cast<size_t>(m_loadStart[0]), event_index,
                  numFillRanges(static_cast<size_t>(m_loadSize[0]), static_cast<size_t>(m_max_id - m_min_id) + 1))) {
    // the processing tasks this replaces would have reported 3 steps each
    prog->reportIncrement(m_loader.splitProcessing ? 6 : 3, entry_name + ": filled events");
    thisBankPulseTimes.reset();
    return;
  }
  if (m_loader.alg->getCancel())
    return;
  // only split if told to and the section to load is at least 1/4 the size
  // of the whole bank
  const bool splitBank = m_loader.splitProcessing && m_max_id > (m_min_id + (bank_size / 4));

  // count the events on each detector ID once, in parallel. This chooses where to split the bank and lets the
  // processing tasks reserve memory without each scanning all the events again.
  std::shared_ptr<std::vector<size_t> const> eventsPerDetId;
  if (splitBank || (m_loader.precount && !useCompressed)) {
    eventsPerDetId = std::make_shared<std::vector<size_t> const>(
        countEventsPerDetId(*event_id, static_cast<size_t>(m_loadSize[0]), m_min_id, m_max_id));
  }

  // split where each half has about the same number of events, so both processing tasks finish together
  const auto mid_id = splitBank ? balancedMidId(*eventsPerDetId, m_min_id) : m_max_id;

  // No error? Launch a new task to process that data.
  const auto numEvents = static_cast<size_t>(m_loadSize[0]);
  const auto startAt = static_cast<size_t>(m_loadStart[0]);

  if (useCompressed) {
    // this method is for unweighted events that the user wants compressed on load

    // TODO should this be created elsewhere?
    const auto [tof_min, tof_max] = Mantid::Kernel::parallel_minmax(event_time_of_flight);

    const bool log_compression = (m_loader.alg->compressTolerance < 0);

    // reduce tof range if filtering was requested
    auto tof_min_fixed = tof_min;
    auto tof_max_fixed = tof_max;
    if (m_loader.alg->filter_tof_range) {
      if (m_loader.alg->filter_tof_max != EMPTY_DBL())
        tof_max_fixed = std::min<float>(tof_max_fixed, static_cast<float>(m_loader.alg->filter_tof_max));
      if (m_loader.alg->filter_tof_min != EMPTY_DBL())
        tof_min_fixed = std::max<float>(tof_min_fixed, static_cast<float>(m_loader.alg->filter_tof_min));
    }

    // fixup the minimum tof for log binning since it cannot be <= 0
    if (log_compression && tof_min_fixed <= 0.) {
      tof_min_fixed = std::abs(static_cast<float>(m_loader.alg->compressTolerance));
    }

    // Join back up the tof limits to the global ones
    // TODO count the bad times-of-flight and discarded events too
    TofStats tofStats;
    tofStats.shortestTof = tof_min_fixed;
    tofStats.longestTof = tof_max_fixed;
    m_loader.addTofStats(tofStats);

    // delta >= 0 is linear, < 0 is log
    double delta = m_loader.alg->compressTolerance;

    // make a vector of logorithmic bins
    auto histogram_bin_edges = std::make_shared<std::vector<double>>();
    std::vector<double> const params{tof_min_fixed, delta, tof_max_fixed + std::abs(delta)};
    Mantid::Kernel::VectorHelper::createAxisFromRebinParams(params, *histogram_bin_edges);

    // create the tasks
    std::shared_ptr<Task> newTask1 = std::make_shared<ProcessBankCompressed>(
        m_loader, entry_name, prog, event_id, event_time_of_flight, startAt, event_index, thisBankPulseTimes, m_min_id,
        mid_id, histogram_bin_edges, m_loader.alg->compressTolerance);
    scheduler.push(newTask1);
    if (m_loader.splitProcessing && (mid_id < m_max_id)) {
      std::shared_ptr<Task> newTask2 = std::make_shared<ProcessBankCompressed>(
          m_loader, entry_name, prog, event_id, event_time_of_flight, startAt, event_index, thisBankPulseTimes,
          (mid_id + 1), m_max_id, histogram_bin_edges, m_loader.alg->compressTolerance);
      scheduler.push(newTask2);
    }
  } else {
    // create all events using traditional method
    std::shared_ptr<Task> newTask1 = std::make_shared<ProcessBankData>(
        m_loader, entry_name, prog, event_id, event_time_of_flight, numEvents, startAt, event_index, thisBankPulseTimes,
        m_have_weight, event_weight, m_min_id, mid_id, eventsPerDetId, m_min_id);
    scheduler.push(newTask1);
    if (m_loader.splitProcessing && (mid_id < m_max_id)) {
      std::shared_ptr<Task> newTask2 = std::make_shared<ProcessBankData>(
          m_loader, entry_name, prog, event_id, event_time_of_flight, numEvents, startAt, event_index,
          thisBankPulseTimes, m_have_weight, event_weight, (mid_id + 1), m_max_id, eventsPerDetId, m_min_id);
      scheduler.push(newTask2);
    }
  }

#ifndef _WIN32
  if (m_loader.alg->getLogger().isDebug())
    m_loader.alg->getLogger().debug() << "Time to LoadBankFromDisk " << entry_name << " " << timer << "\n";
#endif
  thisBankPulseTimes.reset();
}

/** Write the bank's events in [m_min_id, m_max_id] straight into their event lists, in numRanges contiguous ranges of
 * events filled in parallel.
 *
 * A counting pass first finds how many of each range's events go to each list, applying the same pulse and
 * time-of-flight filters as the fill. Each list is then grown once to fit, and each range writes from its own start in
 * every list, right after the earlier ranges' events. The lists therefore get their events in the same order as
 * appending them one at a time would give, with no gaps.
 *
 * This does nothing, and returns false, when two detector IDs with events share one list, as each detector ID gets its
 * own write position, or when the algorithm is cancelled.
 *
 * @param ids :: the detector IDs, where element 0 is event firstEvent of the bank
 * @param tofs :: the times-of-flight in microseconds, matching ids
 * @param firstEvent :: the bank's index of element 0 of the arrays
 * @param event_index :: the bank's event_index
 * @param numRanges :: how many ranges of events to fill in parallel
 * @return true if the events were written
 */
bool LoadBankFromDiskTask::fillInPlace(std::vector<uint32_t> const &ids, std::vector<float> const &tofs,
                                       const size_t firstEvent,
                                       const std::shared_ptr<std::vector<uint64_t>> &event_index,
                                       const size_t numRanges) {
  const auto &lists = m_loader.eventVectors[0];
  const uint32_t minId = m_min_id;
  const uint32_t maxId = m_max_id;
  const auto width = static_cast<size_t>(maxId - minId) + 1;
  const auto tofFilter = m_loader.tofFilter();
  const auto isKept = [&](const uint32_t id, const double tof) {
    return id >= minId && id <= maxId && tofFilter.keeps(tof);
  };
  const auto &pulseTimes = *thisBankPulseTimes;

  const PulseIndexer indexer(event_index, firstEvent, ids.size(), entry_name, m_loader.pulseIndicesToLoad(pulseTimes));
  const PulseRanges ranges(indexer, numRanges);
  const auto rangeCounts = countEventsPerRange(ranges, pulseTimes, ids, tofs, isKept, lists, minId, width);
  std::vector<size_t> totals(width, 0);
  for (const auto &counts : rangeCounts)
    std::transform(totals.cbegin(), totals.cend(), counts.cbegin(), totals.begin(), std::plus<size_t>());
  // each detector ID gets its own write positions, which would collide in a shared list
  if (listsAreShared(lists, minId, totals) || m_loader.alg->getCancel())
    return false;

  auto rangeCursors = growLists(lists, minId, totals, rangeCounts);
  m_loader.addTofStats(fillRanges(ranges, pulseTimes, ids, tofs, isKept, minId, rangeCursors));

  // set the sort order of the lists written, as ProcessBankData does
  const auto sortOrder = pulseTimes.arePulseTimesIncreasing() ? DataObjects::PULSETIME_SORT : DataObjects::UNSORTED;
  auto &outputWS = m_loader.m_ws;
  const size_t numEventLists = outputWS.getNumberHistograms();
  for (size_t i = 0; i < width; ++i) {
    if (totals[i] == 0)
      continue;
    const auto wi = m_loader.workspaceIndexOf(static_cast<detid_t>(minId + i));
    if (wi && *wi < numEventLists)
      outputWS.getSpectrum(*wi).setSortOrder(sortOrder);
  }
  return true;
}

/**
 * Interpret the value describing the number of events. If the number is
 * positive return it unchanged.
 * If the value is negative (can happen at ISIS) add 2^32 to it.
 * @param size :: The size of events value.
 */
uint64_t LoadBankFromDiskTask::recalculateDataSize(const int64_t size) {
  uint64_t ret(size);
  if (size < 0) {
    uint64_t const shift = uint64_t(1) << 32;
    ret += shift;
  }
  return ret;
}

} // namespace Mantid::DataHandling
