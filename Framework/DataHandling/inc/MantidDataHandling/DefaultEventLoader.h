// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/Axis.h"
#include "MantidDataHandling/DllConfig.h"
#include "MantidDataHandling/EventWorkspaceCollection.h"

#include <algorithm>
#include <limits>
#include <mutex>
#include <optional>
#include <string>

class BankPulseTimes;

namespace Mantid {
namespace DataHandling {
class BankPulseTimes;
class LoadEventNexus;

/// The time-of-flight range the user asked to load
class TofFilter {
public:
  TofFilter(const bool filtering, const double tofMin, const double tofMax)
      : m_filtering(filtering), m_tofMin(tofMin), m_tofMax(tofMax) {}
  /// Whether an event with this time-of-flight is loaded
  bool keeps(const double tof) const { return !m_filtering || (tof - m_tofMin) * (tof - m_tofMax) <= 0.; }

private:
  bool m_filtering;
  double m_tofMin;
  double m_tofMax;
};

/// Time-of-flight limits and counts of the events loaded, gathered by each task and then added to the algorithm's
struct TofStats {
  double shortestTof{static_cast<double>(std::numeric_limits<uint32_t>::max()) * 0.1};
  double longestTof{0.};
  /// Events with a time-of-flight too large to be real
  size_t badTofs{0};
  /// Events on detector IDs that have no event list
  size_t discardedEvents{0};

  /// Record the time-of-flight of an event that passed the filters
  void add(const double tof) {
    // Skip any events that are the cause of bad DAS data (e.g. a negative number in uint32 -> 2.4 billion * 100
    // nanosec = 2.4e8 microsec)
    if (tof < 2e8) {
      if (tof > longestTof)
        longestTof = tof;
      if (tof < shortestTof)
        shortestTof = tof;
    } else {
      ++badTofs;
    }
  }

  void merge(const TofStats &other) {
    shortestTof = std::min(shortestTof, other.shortestTof);
    longestTof = std::max(longestTof, other.longestTof);
    badTofs += other.badTofs;
    discardedEvents += other.discardedEvents;
  }
};

/** Helper class for LoadEventNexus that is specific to the current default
  loading code for NXevent_data entries in Nexus files, in particular
  LoadBankFromDiskTask and ProcessBankData.
*/
class MANTID_DATAHANDLING_DLL DefaultEventLoader {
public:
  static void load(LoadEventNexus *alg, EventWorkspaceCollection &ws, bool haveWeights, bool event_id_is_spec,
                   std::vector<std::string> bankNames, const std::vector<int> &periodLog, const std::string &classType,
                   std::vector<std::size_t> bankNumEvents, const bool oldNeXusFileNames, const bool precount,
                   const int chunk, const int totalChunks);

  /// Flag for dealing with a simulated file
  bool m_haveWeights;

  /// True if the event_id is spectrum no not pixel ID
  bool event_id_is_spec;

  /// whether or not to launch multiple ProcessBankData jobs per bank
  bool splitProcessing;

  /// Do we pre-count the # of events in each pixel ID?
  bool precount;

  /// Offset in the pixelID_to_wi_vector to use.
  detid_t pixelID_to_wi_offset;

  /// Maximum (inclusive) event ID possible for this instrument
  int32_t eventid_max{0};

  /// chunk number
  int chunk;
  /// number of chunks
  int totalChunks;
  /// for multiple chunks per bank
  int firstChunkForBank;
  /// number of chunks per bank
  size_t eventsPerChunk;

  LoadEventNexus *alg;
  EventWorkspaceCollection &m_ws;

  /// Vector where index = event_id; value = ptr to std::vector<TofEvent> in the
  /// event list.
  std::vector<std::vector<std::vector<Mantid::Types::Event::TofEvent> *>> eventVectors;

  /// Vector where index = event_id; value = ptr to std::vector<WeightedEvent>
  /// in the event list.
  std::vector<std::vector<std::vector<Mantid::DataObjects::WeightedEvent> *>> weightedEventVectors;

  /// Vector where index = event_id; value = ptr to std::vector<WeightedEventNoTime> in the
  /// event list.
  std::vector<std::vector<std::vector<Mantid::DataObjects::WeightedEventNoTime> *>> weightedNoTimeEventVectors;

  /// Vector where (index = pixel ID+pixelID_to_wi_offset), value = workspace
  /// index)
  std::vector<size_t> pixelID_to_wi_vector;

  /// One entry of pulse times for each preprocessor
  std::vector<std::shared_ptr<BankPulseTimes>> m_bankPulseTimes;
  /// Guards m_bankPulseTimes, which the bank tasks search and add to concurrently
  std::mutex m_bankPulseTimesMutex;

  /// The pulses to load from a bank, after the time and bad pulse filters; empty to load them all
  std::vector<size_t> pulseIndicesToLoad(const BankPulseTimes &pulseTimes) const;
  /// The time-of-flight range to load
  TofFilter tofFilter() const;
  /// Add a task's time-of-flight limits and counts to the algorithm's. Safe to call from several tasks at once.
  void addTofStats(const TofStats &stats) const;
  /// The workspace index for a detector ID (spectrum number if event_id_is_spec), or nothing if the ID is outside the
  /// map. The index is past the last spectrum when the ID has none.
  std::optional<size_t> workspaceIndexOf(const detid_t id) const;

private:
  DefaultEventLoader(LoadEventNexus *alg, EventWorkspaceCollection &ws, bool haveWeights, bool event_id_is_spec,
                     const size_t numBanks, const bool precount, const int chunk, const int totalChunks);
  std::pair<size_t, size_t> setupChunking(std::vector<std::string> &bankNames, std::vector<std::size_t> &bankNumEvents);
  /// Map detector IDs to event lists.
  template <class T> void makeMapToEventLists(std::vector<std::vector<T>> &vectors);
};

/** Generate a look-up table where the index = the pixel ID of an event
 * and the value = a pointer to the EventList in the workspace
 * @param vectors :: the array to create the map on
 */
template <class T> void DefaultEventLoader::makeMapToEventLists(std::vector<std::vector<T>> &vectors) {
  vectors.resize(m_ws.nPeriods());
  if (event_id_is_spec) {
    // Find max spectrum no
    const auto *ax1 = m_ws.getAxis(1);
    specnum_t maxSpecNo = -std::numeric_limits<specnum_t>::max(); // So that any number will be
                                                                  // greater than this
    for (size_t i = 0; i < ax1->length(); i++) {
      specnum_t spec = ax1->spectraNo(i);
      if (spec > maxSpecNo)
        maxSpecNo = spec;
    }

    // These are used by the bank loader to figure out where to put the events
    // The index of eventVectors is a spectrum number so it is simply resized to
    // the maximum
    // possible spectrum number
    eventid_max = maxSpecNo;
    for (size_t i = 0; i < vectors.size(); ++i) {
      vectors[i].resize(maxSpecNo + 1, nullptr);
    }
    for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
      for (size_t i = 0; i < m_ws.getNumberHistograms(); ++i) {
        const auto &spec = m_ws.getSpectrum(i);
        getEventsFrom(m_ws.getSpectrum(i, period), vectors[period][spec.getSpectrumNo()]);
      }
    }
  } else {
    // To avoid going out of range in the vector, this is the MAX INDEX that can
    // go into it
    eventid_max = static_cast<int32_t>(pixelID_to_wi_vector.size() - 1 - pixelID_to_wi_offset);

    // Make an array where index = pixel ID
    // Set the value to NULL by default
    for (size_t i = 0; i < vectors.size(); ++i) {
      vectors[i].resize(eventid_max + 1, nullptr);
    }

    for (size_t j = 0; j < pixelID_to_wi_vector.size(); j++) {
      size_t wi = pixelID_to_wi_vector[j];
      // Save a POINTER to the vector
      if (wi < m_ws.getNumberHistograms()) {
        for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
          getEventsFrom(m_ws.getSpectrum(wi, period), vectors[period][j - pixelID_to_wi_offset]);
        }
      }
    }
  }
}

} // namespace DataHandling
} // namespace Mantid
