#ifndef CHUNKED_JACOBIAN_H
#define CHUNKED_JACOBIAN_H

#include <algorithm>
#include <vector>

#include "autodiff.h"
#include "vector.h"
#include "matrix.h"

namespace ASC_ode {

  using namespace ASC_bla;

  // Number of Jacobian columns one forward-mode pass computes.
  // A pass costs roughly (fixed overhead) + (chunk width) * (work per column),
  // and nested AutoDiff values grow with the width. 16 measured best on chains of 1 to 16 bodies.
#ifndef RIGID_JACOBIAN_CHUNK
#define RIGID_JACOBIAN_CHUNK 16
#endif
  constexpr size_t kJacobianChunk = RIGID_JACOBIAN_CHUNK;

  /**
   * @brief Sorted list of the columns a residual can depend on, without duplicates.
   *
   * Seeding the same column twice in one pass would add up two derivative slots,
   * so duplicates have to go.
   */
  class ColumnSet {
    std::vector<size_t> cols_;
  public:
    void AddRange(size_t first, size_t count) {
      for (size_t i = 0; i < count; i++) cols_.push_back(first + i);
    }
    const std::vector<size_t>& Finish() {
      std::sort(cols_.begin(), cols_.end());
      cols_.erase(std::unique(cols_.begin(), cols_.end()), cols_.end());
      return cols_;
    }
  };

  /**
   * @brief Forward-mode Jacobian of residual(x, f), computed C columns per pass.
   *
   * Only the listed columns are seeded; every other column of df is left untouched,
   * so the caller zeroes df first. residual is called as residual(VectorView<AD>, VectorView<AD>).
   */
  template <size_t C = kJacobianChunk, typename F>
  void ChunkedJacobian(VectorView<double> x, const std::vector<size_t>& cols,
                       size_t dimF, F&& residual, MatrixView<double> df)
  {
    using AD = AutoDiff<C, double>;
    Vector<AD> xd(x.Size());
    for (size_t i = 0; i < x.Size(); i++) xd(i) = AD(x(i));
    Vector<AD> fd(dimF);

    for (size_t start = 0; start < cols.size(); start += C) {
      const size_t n = std::min(C, cols.size() - start);
      for (size_t k = 0; k < n; k++) xd(cols[start + k]).DValue(k) = 1;

      residual(VectorView<AD>(xd), VectorView<AD>(fd));

      for (size_t r = 0; r < dimF; r++)
        for (size_t k = 0; k < n; k++) df(r, cols[start + k]) = fd(r).DValue(k);

      for (size_t k = 0; k < n; k++) xd(cols[start + k]).DValue(k) = 0;
    }
  }

}

#endif
