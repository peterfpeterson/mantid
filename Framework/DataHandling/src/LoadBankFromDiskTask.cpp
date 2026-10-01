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
#include "MantidKernel/ParallelMinMax.h"
#include "MantidKernel/Timer.h"
#include "MantidKernel/Unit.h"
#include "MantidKernel/VectorHelper.h"
#include "MantidNexus/NexusException.h"
#include "MantidNexus/NexusFile.h"
#include "MantidNexus/NexusIOHelper.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <numeric>
#include <sstream>
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

/// Number of events on detector IDs [minId, lastId], where index 0 of counts is minId
size_t eventsUpTo(std::vector<size_t> const &counts, const uint32_t minId, const uint32_t lastId) {
  return std::accumulate(counts.cbegin(), counts.cbegin() + (lastId - minId) + 1, static_cast<size_t>(0));
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

  // Now, we look through existing ones to see if it is already loaded
  // thisBankPulseTimes = NULL;
  for (auto &bankPulseTime : m_loader.m_bankPulseTimes) {
    if (bankPulseTime->equals(thispulseTimes, thisStartTime)) {
      thisBankPulseTimes = bankPulseTime;
      return;
    }
  }

  // Not found? Need to load and add it
  thisBankPulseTimes = std::make_shared<BankPulseTimes>(file, m_framePeriodNumbers);
  m_loader.m_bankPulseTimes.emplace_back(thisBankPulseTimes);
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

  // diagnostics: seconds spent in each phase of reading the bank
  m_loader.diagLog("disk-start bank=" + entry_name + " cost=" + std::to_string(static_cast<size_t>(m_cost)));
  auto diagMark = std::chrono::steady_clock::now();
  auto diagLap = [&diagMark]() {
    const auto now = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(now - diagMark).count();
    diagMark = now;
    return seconds;
  };
  double tOpen = 0., tIndex = 0., tPulse = 0., tId = 0., tTof = 0., tWeights = 0., tClose = 0.;

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
    tOpen = diagLap();

    const bool needPulseInfo = (!m_loader.alg->compressEvents) || m_loader.alg->compressTolerance == 0 ||
                               m_loader.m_ws.nPeriods() > 1 || m_loader.alg->m_is_time_filtered ||
                               m_loader.alg->filter_bad_pulses || m_have_weight;

    // Load the event_index field.
    if (needPulseInfo)
      event_index = this->loadEventIndex(file);
    else
      event_index = nullptr;
    tIndex = diagLap();

    if (!m_loadError) {
      // Load and validate the pulse times
      if (needPulseInfo)
        this->loadPulseTimes(file);
      else
        thisBankPulseTimes = nullptr;
      tPulse = diagLap();
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
        tId = diagLap();

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
          tTof = diagLap();
          if (m_have_weight) {
            event_weight = this->loadEventWeights(file);
            tWeights = diagLap();
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
  tClose = diagLap();
  {
    std::ostringstream msg;
    msg << "disk-read-done bank=" << entry_name << " events=" << m_loadSize[0] << " error=" << m_loadError
        << " open=" << tOpen << " index=" << tIndex << " pulse=" << tPulse << " id=" << tId << " tof=" << tTof
        << " weights=" << tWeights << " close=" << tClose
        << " total=" << (tOpen + tIndex + tPulse + tId + tTof + tWeights + tClose);
    m_loader.diagLog(msg.str());
  }

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
  // only split if told to and the section to load is at least 1/4 the size
  // of the whole bank
  const bool splitBank = m_loader.splitProcessing && m_max_id > (m_min_id + (bank_size / 4));

  // count the events on each detector ID once, in parallel. This chooses where to split the bank and lets the
  // processing tasks reserve memory without each scanning all the events again.
  std::shared_ptr<std::vector<size_t> const> eventsPerDetId;
  double countSeconds = 0.;
  if (splitBank || (m_loader.precount && !useCompressed)) {
    const auto countStart = std::chrono::steady_clock::now();
    eventsPerDetId = std::make_shared<std::vector<size_t> const>(
        countEventsPerDetId(*event_id, static_cast<size_t>(m_loadSize[0]), m_min_id, m_max_id));
    countSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - countStart).count();
  }

  auto mid_id = m_max_id;
  if (splitBank) {
    const auto &counts = *eventsPerDetId;
    const auto midpoint_id = (m_max_id + m_min_id) / 2;
    // split where each half has about the same number of events, so both processing tasks finish together
    mid_id = balancedMidId(counts, m_min_id);

    const size_t total = eventsUpTo(counts, m_min_id, m_max_id);
    const size_t balancedFirst = eventsUpTo(counts, m_min_id, mid_id);
    const size_t midpointFirst = eventsUpTo(counts, m_min_id, midpoint_id);
    std::ostringstream msg;
    msg << "balance bank=" << entry_name << " count_seconds=" << countSeconds << " detids=" << m_min_id << "-"
        << m_max_id << " events=" << total << " midpoint_id=" << midpoint_id << " midpoint_halves=" << midpointFirst
        << "/" << (total - midpointFirst) << " balanced_id=" << mid_id << " balanced_halves=" << balancedFirst << "/"
        << (total - balancedFirst);
    m_loader.diagLog(msg.str());
  } else if (eventsPerDetId) {
    m_loader.diagLog("balance bank=" + entry_name + " count_seconds=" + std::to_string(countSeconds) + " not split");
  }

  // No error? Launch a new task to process that data.
  auto diagQueued = [this](const std::string &type, Task &task, const uint32_t minId, const uint32_t maxId) {
    m_loader.diagLog("queued-process type=" + type + " bank=" + entry_name + " detids=" + std::to_string(minId) + "-" +
                     std::to_string(maxId) + " cost=" + std::to_string(static_cast<size_t>(task.cost())));
  };
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
    // This is not thread safe, so only one thread at a time runs this.
    {
      std::lock_guard<std::mutex> _lock(m_loader.alg->m_tofMutex);
      if (tof_min_fixed < m_loader.alg->shortest_tof) {
        m_loader.alg->shortest_tof = tof_min_fixed;
      }
      if (tof_max_fixed > m_loader.alg->longest_tof) {
        m_loader.alg->longest_tof = tof_max_fixed;
      }
      // TODO
      // m_loader.alg->bad_tofs += badTofs;
      // m_loader.alg->discarded_events += my_discarded_events;
    }

    // delta >= 0 is linear, < 0 is log
    double delta = m_loader.alg->compressTolerance;

    // make a vector of logorithmic bins
    auto histogram_bin_edges = std::make_shared<std::vector<double>>();
    std::vector<double> const params{tof_min_fixed, delta, tof_max_fixed + std::abs(delta)};
    Mantid::Kernel::VectorHelper::createAxisFromRebinParams(params, *histogram_bin_edges);

    // create the tasks
    const std::string TASK_TYPE("compressed");
    std::shared_ptr<Task> newTask1 = std::make_shared<ProcessBankCompressed>(
        m_loader, entry_name, prog, event_id, event_time_of_flight, startAt, event_index, thisBankPulseTimes, m_min_id,
        mid_id, histogram_bin_edges, m_loader.alg->compressTolerance);
    diagQueued(TASK_TYPE, *newTask1, m_min_id, mid_id);
    scheduler.push(newTask1);
    if (m_loader.splitProcessing && (mid_id < m_max_id)) {
      std::shared_ptr<Task> newTask2 = std::make_shared<ProcessBankCompressed>(
          m_loader, entry_name, prog, event_id, event_time_of_flight, startAt, event_index, thisBankPulseTimes,
          (mid_id + 1), m_max_id, histogram_bin_edges, m_loader.alg->compressTolerance);
      diagQueued(TASK_TYPE, *newTask2, mid_id + 1, m_max_id);
      scheduler.push(newTask2);
    }
  } else {
    // create all events using traditional method
    const std::string TASK_TYPE("data");
    std::shared_ptr<Task> newTask1 = std::make_shared<ProcessBankData>(
        m_loader, entry_name, prog, event_id, event_time_of_flight, numEvents, startAt, event_index, thisBankPulseTimes,
        m_have_weight, event_weight, m_min_id, mid_id, eventsPerDetId, m_min_id);
    diagQueued(TASK_TYPE, *newTask1, m_min_id, mid_id);
    scheduler.push(newTask1);
    if (m_loader.splitProcessing && (mid_id < m_max_id)) {
      std::shared_ptr<Task> newTask2 = std::make_shared<ProcessBankData>(
          m_loader, entry_name, prog, event_id, event_time_of_flight, numEvents, startAt, event_index,
          thisBankPulseTimes, m_have_weight, event_weight, (mid_id + 1), m_max_id, eventsPerDetId, m_min_id);
      diagQueued(TASK_TYPE, *newTask2, mid_id + 1, m_max_id);
      scheduler.push(newTask2);
    }
  }

#ifndef _WIN32
  if (m_loader.alg->getLogger().isDebug())
    m_loader.alg->getLogger().debug() << "Time to LoadBankFromDisk " << entry_name << " " << timer << "\n";
#endif
  m_loader.diagLog("disk-end bank=" + entry_name);
  thisBankPulseTimes.reset();
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
