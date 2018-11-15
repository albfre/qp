#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <vector>
#include <cassert>
#include <numeric>

namespace MathUtil {
  using Vector = std::vector<double>;
  using Matrix = std::vector<Vector>;

  void solveLinearSystemOfEquations(const Matrix& A, Vector& b);

  inline double dot(const Vector& a, const Vector& b) {
    assert(a.size() == b.size());
    return std::inner_product(a.cbegin(), a.cend(), b.cbegin(), 0.0);
  }

  inline void matrixTimesVector(const Matrix& m, const Vector& v, Vector& r)
  {
    const auto size = m.size();
    assert(r.size() == size);
    for (size_t i = 0; i < size; ++i) {
      r[i] = dot(m[i], v);
    }
  }

  Matrix matrixTranspose(const Matrix& m);

  template<typename Op>
  void vectorEqOpSelfVector(Vector& r, const Vector& v, const Op& op) {
    const auto size = r.size();
    assert(v.size() == size);
    for (size_t i = 0; i < size; ++i) {
      r[i] = op(r[i], v[i]);
    }
  }

  void vectorMinusEqVector(Vector& r, const Vector& v);
  void vectorPlusEqVector(Vector& r, const Vector& v);
  void vectorPlusEqScalarTimesVector(Vector& r, const double scalar, const Vector& v);

  template<class Op>
  void matrixEqOpMatrixWithOffset(Matrix& M,
                                  const Matrix& M2,
                                  const size_t iOffset,
                                  const size_t jOffset,
                                  const Op& op) {
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
        mRow[jOffset + j] = op(m2Row[j]);
      }
    }
  }
}

#endif //MATHUTIL_H
