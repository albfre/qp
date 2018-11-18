#include <cmath>
#include <functional>
#include <iostream>
#include "MathUtil.h"

namespace MathUtil {
  void vectorMinusEqVector(Vector& r, const Vector& v) { vectorEqOpSelfVector(r, v, std::minus<double>()); }
  void vectorPlusEqVector(Vector& r, const Vector& v) { vectorEqOpSelfVector(r, v, std::plus<double>()); }
  void vectorPlusEqScalarTimesVector(Vector& r, const double scalar, const Vector& v) { vectorEqOpSelfVector(r, v, [scalar] (const auto& rx, const auto& vx) { return rx + scalar * vx; }); }

  Matrix matrixTranspose(const Matrix& m)
  {
    const auto ySize = m.size();
    assert(ySize > 0);
    const auto xSize = m.front().size();
    Matrix mT(xSize, Vector(ySize));
    for (size_t i = 0; i < ySize; ++i) {
      const auto& row = m[i];
      assert(row.size() == xSize);
      for (size_t j = 0; j < xSize; ++j) {
        mT[j][i] = row[j];
      }
    }
    return mT;
  }

  void matrixEqMatrixWithOffset(Matrix& M,
                                const Matrix& M2,
                                const size_t iOffset,
                                const size_t jOffset) {
    if (M2.empty()) {
      return;
    }
    const auto m = M2.size();
    assert(M.size() >= iOffset + m);
    const auto n = M2.front().size();
    assert(M.front().size() >= jOffset + n);

    for (size_t i = 0; i < m; ++i) {
      const auto& m2Row = M2[i];
      assert(m2Row.size() == n);
      auto& mRow = M[iOffset + i];
      assert(mRow.size() >= jOffset + n);
      for (size_t j = 0; j < n; ++j) {
        mRow[jOffset + j] = m2Row[j];
      }
    }
  }
}
