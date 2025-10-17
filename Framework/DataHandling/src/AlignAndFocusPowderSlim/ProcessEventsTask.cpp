// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "ProcessEventsTask.h"
#include <span>

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

ProcessEventsTask::ProcessEventsTask(const std::vector<float> *tofs, const std::vector<double> *binedges)
    : y_temp(binedges->size() - 1, 0), m_tofs(tofs), m_binedges(binedges) {}

ProcessEventsTask::ProcessEventsTask(ProcessEventsTask &other, tbb::split)
    : y_temp(other.y_temp.size(), 0), m_tofs(other.m_tofs), m_binedges(other.m_binedges) {}

void ProcessEventsTask::operator()(const tbb::blocked_range<size_t> &range) {
  // Cache values to reduce number of function calls
  const auto &binedges_cbegin = m_binedges->cbegin();
  const auto &binedges_cend = m_binedges->cend();
  const auto &tof_min = m_binedges->front();

  std::span<const float> tofs(m_tofs->begin() + range.begin(), range.size());
  for (const auto tof : tofs) {
    const double tof_d = static_cast<double>(tof);
    if (!(tof_d < tof_min)) { // check against max done in transformation
      // Find the bin index using binary search
      const auto &it = std::upper_bound(binedges_cbegin, binedges_cend, tof_d);

      // Increment the count if a bin was found
      const auto &bin = static_cast<size_t>(std::distance(binedges_cbegin, it) - 1);
      y_temp[bin]++;
    }
  }
}

void ProcessEventsTask::join(const ProcessEventsTask &other) {
  // Combine local histograms
  std::transform(y_temp.begin(), y_temp.end(), other.y_temp.cbegin(), y_temp.begin(), std::plus<>{});
}

} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
