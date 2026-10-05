// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include <utility>

#include "MantidDataHandling/DefaultEventLoader.h"
#include "MantidDataHandling/LoadEventNexus.h"
#include "MantidDataHandling/ProcessBankData.h"
#include "MantidDataHandling/PulseIndexer.h"
#include "MantidKernel/Timer.h"

#include <chrono>
#include <sstream>
#ifdef __linux__
#include <sys/resource.h>
#endif

using namespace Mantid::DataObjects;

namespace Mantid::DataHandling {

namespace {
/// Diagnostics: minor page faults of the calling thread so far, or 0 where that is not available
long threadMinorFaults() {
#ifdef __linux__
  rusage usage{};
  if (getrusage(RUSAGE_THREAD, &usage) == 0)
    return usage.ru_minflt;
#endif
  return 0;
}

/// Diagnostics: wall time and minor page faults of the calling thread between laps
class DiagLap {
public:
  DiagLap() : m_time(std::chrono::steady_clock::now()), m_faults(threadMinorFaults()) {}
  /// Seconds and minor page faults since the previous lap, or since construction
  std::pair<double, long> lap() {
    const auto now = std::chrono::steady_clock::now();
    const auto faults = threadMinorFaults();
    const std::pair<double, long> result{std::chrono::duration<double>(now - m_time).count(), faults - m_faults};
    m_time = now;
    m_faults = faults;
    return result;
  }

private:
  std::chrono::steady_clock::time_point m_time;
  long m_faults;
};
} // namespace

ProcessBankData::ProcessBankData(DefaultEventLoader &m_loader, const std::string &entry_name, API::Progress *prog,
                                 std::shared_ptr<std::vector<uint32_t>> const &event_id,
                                 std::shared_ptr<std::vector<float>> const &tevent_time_of_flight, size_t numEvents,
                                 size_t startAt, std::shared_ptr<std::vector<uint64_t>> const &tevent_index,
                                 std::shared_ptr<BankPulseTimes> const &thisBankPulseTimes, bool have_weight,
                                 std::shared_ptr<std::vector<float>> const &tevent_weight, detid_t min_event_id,
                                 detid_t max_event_id, std::shared_ptr<std::vector<size_t> const> eventsPerDetId,
                                 detid_t eventsPerDetIdMin)
    : Task(), m_loader(m_loader), entry_name(std::move(entry_name)), prog(prog), event_detid(event_id),
      event_time_of_flight(tevent_time_of_flight), numEvents(numEvents), startAt(startAt), event_index(tevent_index),
      thisBankPulseTimes(thisBankPulseTimes), have_weight(have_weight), event_weight(tevent_weight),
      m_min_detid(min_event_id), m_max_detid(max_event_id), m_eventsPerDetId(std::move(eventsPerDetId)),
      m_eventsPerDetIdMin(eventsPerDetIdMin) {
  // Cost is approximately proportional to the number of events to process.
  m_cost = static_cast<double>(numEvents);

  if (m_max_detid < m_min_detid) {
    std::stringstream msg;
    msg << "max detid (" << m_max_detid << ") < min (" << m_min_detid << ")";
    throw std::runtime_error(msg.str());
  }
}

/*
 * Pre-counting the events per pixel ID allows for allocating the proper amount of memory in each output event vector
 */
void ProcessBankData::preCountAndReserveMem() {
  // ---- Pre-counting events per pixel ID ----
  DiagLap diagLap;
  // use the counts shared by this bank's tasks when they cover this range, rather than scanning all the events again
  if (!eventsPerDetIdInRange()) {
    m_localCounts.assign(m_max_detid - m_min_detid + 1, 0);
    for (size_t i = 0; i < numEvents; i++) {
      const auto thisId = static_cast<detid_t>((*event_detid)[i]);
      if (!(thisId < m_min_detid || thisId > m_max_detid)) // or allows for skipping out early
        m_localCounts[thisId - m_min_detid]++;
    }
  }
  // index 0 is m_min_detid
  const size_t *counts = eventsPerDetIdInRange();
  std::tie(m_diagCountSeconds, m_diagCountFaults) = diagLap.lap();

  // Now we pre-allocate (reserve) the vectors of events in each pixel counted
  auto &outputWS = m_loader.m_ws;
  const auto *alg = m_loader.alg;
  const size_t numEventLists = outputWS.getNumberHistograms();
  for (detid_t pixID = m_min_detid; pixID <= m_max_detid; ++pixID) {
    const auto pixelIndex = pixID - m_min_detid; // index from zero
    if (counts[pixelIndex] > 0) {
      const size_t wi = getWorkspaceIndexFromPixelID(pixID);
      // Find the workspace index corresponding to that pixel ID
      // Allocate it
      if (wi < numEventLists) {
        outputWS.reserveEventListAt(wi, counts[pixelIndex]);
      }
      if ((wi % 20 == 0) && alg->getCancel())
        return; // User cancellation
    }
  }
  std::tie(m_diagReserveSeconds, m_diagReserveFaults) = diagLap.lap();
}

/** Run the data processing
 * FIXME/TODO - split run() into readable methods
 */
void ProcessBankData::run() {
  // timer for performance
  Mantid::Kernel::Timer timer;
  const std::string diagName =
      "type=data bank=" + entry_name + " detids=" + std::to_string(m_min_detid) + "-" + std::to_string(m_max_detid);
  m_loader.diagLog("process-start " + diagName + " cost=" + std::to_string(static_cast<size_t>(m_cost)));
  const double diagStart = m_loader.diagElapsed();
  // diagnostics: per-phase wall time and minor page faults, and how many events and detector IDs this task fills
  std::pair<double, long> diagSetup, diagFill, diagFinal, diagTofJoin;
  size_t diagEventsAdded = 0;
  size_t diagDetIdsTouched = 0;

  // Local tof limits and counts
  TofStats tofStats;

  prog->report(entry_name + ": precount");
  // ---- Pre-counting events per pixel ID ----
  if (m_loader.precount) {
    this->preCountAndReserveMem();
    if (m_loader.alg->getCancel())
      return; // User cancellation
  }

  DiagLap diagLap;
  // this assumes that pulse indices are sorted
  if (!std::is_sorted(event_index->cbegin(), event_index->cend()))
    throw std::runtime_error("Event index is not sorted");

  // And there are this many pulses
  prog->report(entry_name + ": filling events");

  auto *alg = m_loader.alg;

  // Will we need to compress?
  const bool compress = (alg->compressEvents);

  // Which detector IDs were touched? The pre-count already knows which detector IDs have events, so only track them
  // here when there is no pre-count
  const size_t *counts = eventsPerDetIdInRange();
  const bool trackUsedDetIds = (counts == nullptr);
  std::vector<bool> usedDetIds(trackUsedDetIds ? m_max_detid - m_min_detid + 1 : 0, false);

  const auto tofFilter = m_loader.tofFilter();

  // set up wall-clock filtering if it was requested
  const PulseIndexer pulseIndexer(event_index, startAt, numEvents, entry_name,
                                  m_loader.pulseIndicesToLoad(*thisBankPulseTimes));
  diagSetup = diagLap.lap();

  // loop over all pulses
  for (const auto &pulseIter : pulseIndexer) {
    // Save the pulse time at this index for creating those events
    const auto &pulsetime = thisBankPulseTimes->pulseTime(pulseIter.pulseIndex);
    const int logPeriodNumber = thisBankPulseTimes->periodNumber(pulseIter.pulseIndex);
    const auto periodIndex = static_cast<size_t>(logPeriodNumber - 1);

    // loop through events associated with a single pulse
    for (std::size_t eventIndex = pulseIter.eventIndexStart; eventIndex < pulseIter.eventIndexStop; ++eventIndex) {
      // We cached a pointer to the vector<tofEvent> -> so retrieve it and add
      // the event
      const detid_t &detId = static_cast<detid_t>((*event_detid)[eventIndex]);
      if (detId >= m_min_detid && detId <= m_max_detid) {
        // Create the tofevent
        const auto tof = static_cast<double>((*event_time_of_flight)[eventIndex]);
        if (tofFilter.keeps(tof)) {
          // Handle simulated data if present
          if (have_weight) {
            auto *eventVector = m_loader.weightedEventVectors[periodIndex][detId];
            // NULL eventVector indicates a bad spectrum lookup
            if (eventVector) {
              const auto weight = static_cast<double>((*event_weight)[eventIndex]);
              const double errorSq = weight * weight;
              eventVector->emplace_back(tof, pulsetime, weight, errorSq);
              ++diagEventsAdded;
            } else {
              ++tofStats.discardedEvents;
            }
          } else {
            // We have cached the vector of events for this detector ID
            auto *eventVector = m_loader.eventVectors[periodIndex][detId];
            // NULL eventVector indicates a bad spectrum lookup
            if (eventVector) {
              eventVector->emplace_back(std::move(tof), pulsetime);
              ++diagEventsAdded;
            } else {
              ++tofStats.discardedEvents;
            }
          }

          // tof limits from things observed here
          tofStats.add(tof);

          // Track all the touched wi
          if (trackUsedDetIds) {
            const auto detidIndex = detId - m_min_detid;
            if (!usedDetIds[detidIndex])
              usedDetIds[detidIndex] = true;
          }
        } // valid time-of-flight

      } // valid detector IDs
    } // for events in pulse
    // check if cancelled after each 100s of pulses (assumes 60Hz)
    if ((pulseIter.pulseIndex % 6000 == 0) && alg->getCancel())
      return;
  } // for pulses

  diagFill = diagLap.lap();

  // Default pulse time (if none are found)
  const auto pulseSortingType =
      thisBankPulseTimes->arePulseTimesIncreasing() ? DataObjects::PULSETIME_SORT : DataObjects::UNSORTED;

  //------------ Compress Events (or set sort order) ------------------
  // Do it on all the detector IDs we touched
  auto &outputWS = m_loader.m_ws;
  const size_t numEventLists = outputWS.getNumberHistograms();
  for (detid_t pixID = m_min_detid; pixID <= m_max_detid; ++pixID) {
    const auto detidIndex = static_cast<size_t>(pixID - m_min_detid);
    if (trackUsedDetIds ? usedDetIds[detidIndex] : counts[detidIndex] > 0) {
      ++diagDetIdsTouched;
      // Find the workspace index corresponding to that pixel ID
      size_t wi = getWorkspaceIndexFromPixelID(pixID);
      if (wi < numEventLists) {
        auto &el = outputWS.getSpectrum(wi);
        // set the sort order based on what is known
        el.setSortOrder(pulseSortingType);
        // compress events if requested
        if (compress)
          el.compressEvents(alg->compressTolerance, &el);
      }
    }
  }
  diagFinal = diagLap.lap();
  prog->report(entry_name + ": filled events");

  alg->getLogger().debug() << entry_name << (thisBankPulseTimes->arePulseTimesIncreasing() ? " had " : " DID NOT have ")
                           << "monotonically increasing pulse times\n";

  // Join back up the tof limits to the global ones
  m_loader.addTofStats(tofStats);
  diagTofJoin = diagLap.lap();

#ifndef _WIN32
  if (alg->getLogger().isDebug())
    alg->getLogger().debug() << "Time to ProcessBankData " << entry_name << " " << timer << "\n";
#endif
  event_detid.reset();
  event_time_of_flight.reset();
  event_index.reset();
  event_weight.reset();
  thisBankPulseTimes.reset();
  m_loader.diagLog("process-end " + diagName + " duration=" + std::to_string(m_loader.diagElapsed() - diagStart));
  {
    std::ostringstream msg;
    msg << "process-phases " << diagName << " events_added=" << diagEventsAdded
        << " detids_touched=" << diagDetIdsTouched << " precount_count=" << m_diagCountSeconds << "s/"
        << m_diagCountFaults << "pf precount_reserve=" << m_diagReserveSeconds << "s/" << m_diagReserveFaults
        << "pf setup=" << diagSetup.first << "s/" << diagSetup.second << "pf fill=" << diagFill.first << "s/"
        << diagFill.second << "pf final=" << diagFinal.first << "s/" << diagFinal.second
        << "pf tofjoin=" << diagTofJoin.first << "s/" << diagTofJoin.second << "pf";
    m_loader.diagLog(msg.str());
  }
} // END-OF-RUN()

/**
 * The number of events on each detector ID in [m_min_detid, m_max_detid], from the counts shared by this bank's
 * tasks when they cover this range, or else from this task's own pre-count.
 * @return pointer where index 0 is m_min_detid, or nullptr if the events have not been counted
 */
const size_t *ProcessBankData::eventsPerDetIdInRange() const {
  if (m_eventsPerDetId && m_min_detid >= m_eventsPerDetIdMin &&
      static_cast<size_t>(m_max_detid - m_eventsPerDetIdMin) < m_eventsPerDetId->size())
    return m_eventsPerDetId->data() + (m_min_detid - m_eventsPerDetIdMin);
  if (!m_localCounts.empty())
    return m_localCounts.data();
  return nullptr;
}

/**
 * Get the workspace index for a given pixel ID. Throws if the pixel ID is
 * not in the expected range.
 *
 * @param pixID :: The pixel ID to look up
 * @return The workspace index for this pixel
 */
size_t ProcessBankData::getWorkspaceIndexFromPixelID(const detid_t pixID) {
  const auto wi = m_loader.workspaceIndexOf(pixID);
  if (!wi) {
    std::stringstream msg;
    msg << "Error finding workspace index; pixelID " << pixID << " with offset " << m_loader.pixelID_to_wi_offset
        << " is out of range (length=" << m_loader.pixelID_to_wi_vector.size() << ")";
    throw std::runtime_error(msg.str());
  }
  return *wi;
}
} // namespace Mantid::DataHandling
