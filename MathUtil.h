#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <cassert>
#include <numeric>
#include <vector>

namespace MathUtil {
  using Vector = std::vector<double>;
  using Matrix = std::vector<Vector>;

  inline double dot(const Vector& a, const Vector& b) {
    assert(a.size() == b.size());
    return std::inner_product(a.cbegin(), a.cend(), b.cbegin(), 0.0);
  }

  inline Vector matrixTimesVector(const Matrix& m, const Vector& v)
  {
    const auto size = m.size();
    Vector r(size);
    for (size_t i = 0; i < size; ++i) {
      r[i] = dot(m[i], v);
    }
    return r;
  }

  template<class Op>
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
  void vectorTimesEqVector(Vector& r, const Vector& v);
  void vectorDivEqMinusVector(Vector& r, const Vector& v);

  Matrix matrixTranspose(const Matrix& m);

  void matrixEqMatrixWithOffset(Matrix& M,
                                const Matrix& M2,
                                const size_t iOffset,
                                const size_t jOffset);
}

#endif //MATHUTIL_H
