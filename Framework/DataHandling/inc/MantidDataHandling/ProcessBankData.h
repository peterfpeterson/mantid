// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataHandling/BankPulseTimes.h"
#include "MantidGeometry/IDTypes.h"
#include "MantidKernel/Task.h"

#include <memory>

namespace Mantid {
namespace API {
class Progress;
}
namespace DataHandling {
class DefaultEventLoader;

/** This task does the disk IO from loading the NXS file,
 * and so will be on a disk IO mutex */
class ProcessBankData : public Mantid::Kernel::Task {
public:
  /** Constructor
   *
   * @param loader :: DefaultEventLoader
   * @param entry_name :: name of the bank
   * @param prog :: Progress reporter
   * @param event_id :: array with event IDs
   * @param event_time_of_flight :: array with event TOFS
   * @param numEvents :: how many events in the arrays
   * @param startAt :: index of the first event from event_index
   * @param event_index :: vector of event index (length of # of pulses)
   * @param thisBankPulseTimes :: ptr to the pulse times for this particular
   *bank.
   * @param have_weight :: flag for handling simulated files
   * @param event_weight :: array with weights for events
   * @param min_event_id ;: minimum detector ID to load
   * @param max_event_id :: maximum detector ID to load
   * @param eventsPerDetId :: optional number of events on each detector ID, covering at least
   * [min_event_id, max_event_id]. When supplied, pre-counting uses it rather than scanning the events again.
   * @param eventsPerDetIdMin :: detector ID of index 0 in eventsPerDetId
   */
  ProcessBankData(DefaultEventLoader &loader, const std::string &entry_name, API::Progress *prog,
                  std::shared_ptr<std::vector<uint32_t>> const &event_id,
                  std::shared_ptr<std::vector<float>> const &event_time_of_flight, size_t numEvents, size_t startAt,
                  std::shared_ptr<std::vector<uint64_t>> const &event_index,
                  std::shared_ptr<BankPulseTimes> const &thisBankPulseTimes, bool have_weight,
                  std::shared_ptr<std::vector<float>> const &event_weight, detid_t min_event_id, detid_t max_event_id,
                  std::shared_ptr<std::vector<size_t> const> eventsPerDetId = nullptr, detid_t eventsPerDetIdMin = 0);

  void run() override;

private:
  size_t getWorkspaceIndexFromPixelID(const detid_t pixID);
  void preCountAndReserveMem();
  const size_t *eventsPerDetIdInRange() const;

  /// Algorithm being run
  DefaultEventLoader &m_loader;
  /// NXS address to bank
  std::string entry_name;
  /// Progress reporting
  API::Progress *prog;
  /// event pixel ID array
  std::shared_ptr<std::vector<uint32_t> const> event_detid;
  /// event TOF array
  std::shared_ptr<std::vector<float> const> event_time_of_flight;
  /// # of events in arrays
  size_t numEvents;
  /// index of the first event from event_index
  size_t startAt;
  /// vector of event index (length of # of pulses)
  std::shared_ptr<std::vector<uint64_t> const> event_index;
  /// Pulse times for this bank
  std::shared_ptr<BankPulseTimes const> thisBankPulseTimes;
  /// Flag for simulated data
  bool have_weight;
  /// event weights array
  std::shared_ptr<std::vector<float> const> event_weight;
  /// Minimum pixel id (inclusive)
  detid_t m_min_detid;
  /// Maximum pixel id (inclusive)
  detid_t m_max_detid;
  /// Optional number of events on each detector ID, shared by the tasks for one bank
  std::shared_ptr<std::vector<size_t> const> m_eventsPerDetId;
  /// Detector ID of index 0 in m_eventsPerDetId
  detid_t m_eventsPerDetIdMin;
  /// Number of events on each detector ID in [m_min_detid, m_max_detid], when counted by this task
  std::vector<size_t> m_localCounts;
  /// Diagnostics: seconds and minor page faults spent counting and reserving in preCountAndReserveMem
  double m_diagCountSeconds{0.};
  double m_diagReserveSeconds{0.};
  long m_diagCountFaults{0};
  long m_diagReserveFaults{0};
}; // ENDDEF-CLASS ProcessBankData
} // namespace DataHandling
} // namespace Mantid
