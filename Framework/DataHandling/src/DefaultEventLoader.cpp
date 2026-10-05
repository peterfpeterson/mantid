// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataHandling/DefaultEventLoader.h"
#include "MantidAPI/Progress.h"
#include "MantidDataHandling/BankPulseTimes.h"
#include "MantidDataHandling/LoadBankFromDiskTask.h"
#include "MantidDataHandling/LoadEventNexus.h"
#include "MantidKernel/ThreadPool.h"
#include "MantidKernel/ThreadScheduler.h"
#include "MantidKernel/TimeROI.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <thread>
#ifdef __linux__
#include <unistd.h>
#endif

using namespace Mantid::Kernel;

namespace Mantid::DataHandling {

void DefaultEventLoader::load(LoadEventNexus *alg, EventWorkspaceCollection &ws, bool haveWeights,
                              bool event_id_is_spec, std::vector<std::string> bankNames,
                              const std::vector<int> &periodLog, const std::string &classType,
                              std::vector<std::size_t> bankNumEvents, const bool oldNeXusFileNames, const bool precount,
                              const int chunk, const int totalChunks) {
  DefaultEventLoader loader(alg, ws, haveWeights, event_id_is_spec, bankNames.size(), precount, chunk, totalChunks);

  auto bankRange = loader.setupChunking(bankNames, bankNumEvents);

  // Make the thread pool
  auto scheduler = new ThreadSchedulerLargestCost;
  ThreadPool pool(scheduler);
  loader.m_diagStart = std::chrono::steady_clock::now();
  loader.m_diagScheduler = scheduler;
  loader.diagLog("load-begin banks=" + std::to_string(bankRange.second - bankRange.first) +
                 " threads=" + std::to_string(ThreadPool::getNumPhysicalCores()) + " scheduler=LargestCost");

  // set up progress bar for the rest of the (multi-threaded) process
  // only the non-empty banks in this chunk's range get tasks, so only count those
  const auto numBanksToLoad = static_cast<size_t>(std::count_if(bankNumEvents.cbegin() + bankRange.first,
                                                                bankNumEvents.cbegin() + bankRange.second,
                                                                [](const auto numEvents) { return numEvents > 0; }));
  size_t numProg = numBanksToLoad * (1 + 3); // 1 = disktask, 3 = proc task
  if (loader.splitProcessing)
    numProg += numBanksToLoad * 3; // 3 = second proc task
  auto prog = std::make_unique<API::Progress>(loader.alg, 0.3, 1.0, numProg);

  for (size_t i = bankRange.first; i < bankRange.second; i++) {
    if (bankNumEvents[i] > 0) {
      loader.diagLog("queued-disk bank=" + bankNames[i] + " cost=" + std::to_string(bankNumEvents[i]));
      pool.schedule(std::make_shared<LoadBankFromDiskTask>(loader, bankNames[i], classType, bankNumEvents[i],
                                                           oldNeXusFileNames, prog.get(), *scheduler, periodLog));
    } else {
      loader.diagLog("skipped-empty bank=" + bankNames[i]);
    }
  }
  // Start and end all threads
  pool.joinAll();
  loader.diagLog("load-end");
  loader.m_diagScheduler = nullptr;
}

DefaultEventLoader::DefaultEventLoader(LoadEventNexus *alg, EventWorkspaceCollection &ws, bool haveWeights,
                                       bool event_id_is_spec, const size_t numBanks, const bool precount,
                                       const int chunk, const int totalChunks)
    : m_haveWeights(haveWeights), event_id_is_spec(event_id_is_spec), precount(precount), chunk(chunk),
      totalChunks(totalChunks), firstChunkForBank(1), eventsPerChunk(0), alg(alg), m_ws(ws) {
  // This map will be used to find the workspace index
  if (event_id_is_spec)
    pixelID_to_wi_vector = m_ws.getSpectrumToWorkspaceIndexVector(pixelID_to_wi_offset);
  else
    pixelID_to_wi_vector = m_ws.getDetectorIDToWorkspaceIndexVector(pixelID_to_wi_offset, true);

  // Cache a map for speed.
  if (!haveWeights) {
    if (alg->compressEvents && alg->compressTolerance != 0) {
      // Convert to weighted events
      for (size_t i = 0; i < m_ws.getNumberHistograms(); i++) {
        for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
          m_ws.getSpectrum(i, period).switchTo(API::WEIGHTED_NOTIME);
        }
      }
      makeMapToEventLists(weightedNoTimeEventVectors);
    } else {
      makeMapToEventLists(eventVectors);
    }
  } else {
    // Convert to weighted events
    for (size_t i = 0; i < m_ws.getNumberHistograms(); i++) {
      for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
        m_ws.getSpectrum(i, period).switchTo(API::WEIGHTED);
      }
    }
    makeMapToEventLists(weightedEventVectors);
  }

  // split banks up if the number of cores is more than twice the number of
  // banks
  splitProcessing = bool(numBanks * 2 < ThreadPool::getNumPhysicalCores());
}

std::pair<size_t, size_t> DefaultEventLoader::setupChunking(std::vector<std::string> &bankNames,
                                                            std::vector<std::size_t> &bankNumEvents) {
  size_t bank0 = 0;
  size_t bankn = bankNames.size();
  if (chunk != EMPTY_INT()) // We are loading part - work out the bank number range
  {
    const size_t total_events = std::accumulate(bankNumEvents.cbegin(), bankNumEvents.cend(), static_cast<size_t>(0));
    eventsPerChunk = total_events / totalChunks;
    // Sort banks by size
    size_t tmp;
    std::string stmp;
    for (size_t i = 0; i < bankn; i++)
      for (size_t j = 0; j < bankn - 1; j++)
        if (bankNumEvents[j] < bankNumEvents[j + 1]) {
          tmp = bankNumEvents[j];
          bankNumEvents[j] = bankNumEvents[j + 1];
          bankNumEvents[j + 1] = tmp;
          stmp = bankNames[j];
          bankNames[j] = bankNames[j + 1];
          bankNames[j + 1] = stmp;
        }
    const auto bigBanks = std::count_if(bankNumEvents.cbegin(), bankNumEvents.cend(),
                                        [this](const size_t numEvents) { return numEvents > eventsPerChunk; });
    // Each chunk is part of bank or multiple whole banks
    // 0.5 for last chunk of a bank with multiple chunks
    // 0.1 for multiple whole banks not completely filled
    eventsPerChunk +=
        static_cast<size_t>((static_cast<double>(bigBanks) / static_cast<double>(totalChunks) * 0.5 + 0.05) *
                            static_cast<double>(eventsPerChunk));
    double partialChunk = 0.;
    firstChunkForBank = 1;
    for (int chunki = 1; chunki <= chunk; chunki++) {
      if (partialChunk > 1.) {
        partialChunk = 0.;
        firstChunkForBank = chunki;
        bank0 = bankn;
      }
      if (bankNumEvents[bank0] > 1) {
        partialChunk += static_cast<double>(eventsPerChunk) / static_cast<double>(bankNumEvents[bank0]);
      }
      if (chunki < totalChunks)
        bankn = bank0 + 1;
      else
        bankn = bankNames.size();
      if (chunki == firstChunkForBank && partialChunk > 1.0)
        bankn += static_cast<size_t>(partialChunk) - 1;
      if (bankn > bankNames.size())
        bankn = bankNames.size();
    }
    for (size_t i = bank0; i < bankn; i++) {
      size_t start_event = (chunk - firstChunkForBank) * eventsPerChunk;
      size_t stop_event = bankNumEvents[i];
      // Don't change stop_event for the final chunk
      if (start_event + eventsPerChunk < stop_event)
        stop_event = start_event + eventsPerChunk;
      bankNumEvents[i] = stop_event - start_event;
    }
  }
  return {bank0, bankn};
}

/** The pulses to load from a bank, after the time and bad pulse filters
 * @param pulseTimes :: the bank's pulse times
 * @return the indices of the pulses to load, as the PulseIndexer takes them; empty to load them all
 */
std::vector<size_t> DefaultEventLoader::pulseIndicesToLoad(const BankPulseTimes &pulseTimes) const {
  std::vector<size_t> pulseROI;
  if (alg->m_is_time_filtered)
    pulseROI = pulseTimes.getPulseIndices(alg->filter_time_start, alg->filter_time_stop);
  if (alg->filter_bad_pulses)
    pulseROI = Kernel::ROI::calculate_intersection(
        pulseROI, pulseTimes.getPulseIndices(alg->bad_pulses_timeroi->toTimeIntervals()));
  return pulseROI;
}

TofFilter DefaultEventLoader::tofFilter() const {
  return TofFilter(alg->filter_tof_range, alg->filter_tof_min, alg->filter_tof_max);
}

void DefaultEventLoader::addTofStats(const TofStats &stats) const {
  std::lock_guard<std::mutex> lock(alg->m_tofMutex);
  alg->shortest_tof = std::min(alg->shortest_tof, stats.shortestTof);
  alg->longest_tof = std::max(alg->longest_tof, stats.longestTof);
  alg->bad_tofs += stats.badTofs;
  alg->discarded_events += stats.discardedEvents;
}

std::optional<size_t> DefaultEventLoader::workspaceIndexOf(const detid_t id) const {
  const detid_t offsetId = id + pixelID_to_wi_offset;
  if (offsetId < 0 || static_cast<size_t>(offsetId) >= pixelID_to_wi_vector.size())
    return std::nullopt;
  return pixelID_to_wi_vector[static_cast<size_t>(offsetId)];
}

double DefaultEventLoader::diagElapsed() const {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - m_diagStart).count();
}

void DefaultEventLoader::diagLog(const std::string &msg) const {
  std::ostringstream line;
  line << "[LEN-DIAG] t=" << std::fixed << std::setprecision(4) << diagElapsed()
       << " tid=" << std::this_thread::get_id();
  if (m_diagScheduler)
    line << " queue=" << m_diagScheduler->size();
#ifdef __linux__
  // resident set size in MiB, from the second field of /proc/self/statm (in pages)
  {
    std::ifstream statm("/proc/self/statm");
    size_t sizePages = 0, residentPages = 0;
    if (statm >> sizePages >> residentPages)
      line << " rss_mib=" << (residentPages * static_cast<size_t>(sysconf(_SC_PAGESIZE))) / (1024 * 1024);
  }
#endif
  line << " " << msg << "\n";
  alg->getLogger().notice() << line.str();
}

} // namespace Mantid::DataHandling
