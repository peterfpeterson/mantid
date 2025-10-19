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

  // modeled from EventList::findFirstEvent
  auto tof_iter = std::find_if_not(tofs.begin(), tofs.end(),
                                   [tof_min](const float tof) { return tof < static_cast<float>(tof_min); });
  const auto tof_iter_end = tofs.end();

  // go through all the events - modeled from EventList::generateCountsHistogram
  for (auto itx = binedges_cbegin; tof_iter != tof_iter_end; ++tof_iter) {
    const auto tof = static_cast<double>(*tof_iter);
    itx = std::find_if(itx, binedges_cend, [tof](const double edge) { return !(tof < edge); });
    if (itx == binedges_cend) {
      break; // all done
    }
    const auto bin = static_cast<size_t>(std::max(std::distance(binedges_cbegin, itx) - 1, std::ptrdiff_t{0}));
    y_temp[bin]++;
  }

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
