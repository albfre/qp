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

  inline Matrix matrixTimesMatrix(const Matrix& m1, const Matrix& m2)
  {
    const auto m = m1.size();
    const auto n = m2.size();
    assert(m > 0);
    assert(n > 0);
    assert(m1.front().size() == n);
    const auto p = m2.front().size();
    Matrix mat(m, Vector(p));
    for (size_t i = 0; i < m; ++i) {
      const auto& v = m1[i];
      for (size_t j = 0; j < p; ++j) {
        for (size_t k = 0; k < n; ++k) {
          mat[i][j] += v[k] * m2[k][j];
        }
      }
    }

    return mat;
  }

  template<class Op>
  void vectorEqOpSelfVector(Vector& r, const Vector& v, const Op& op) {
    const auto size = r.size();
    assert(v.size() == size);
    for (size_t i = 0; i < size; ++i) {
      r[i] = op(r[i], v[i]);
    }
  }

  void vectorPlusEqScalar(Vector& r, double scalar);
  void vectorTimesEqScalar(Vector& r, double scalar);

  void vectorMinusEqVector(Vector& r, const Vector& v);
  void vectorPlusEqVector(Vector& r, const Vector& v);
  void vectorPlusEqScalarTimesVector(Vector& r, double scalar, const Vector& v);
  void vectorTimesEqVector(Vector& r, const Vector& v);
  void vectorDivEqMinusVector(Vector& r, const Vector& v);

  Matrix matrixTranspose(const Matrix& m);

  void matrixEqMatrixWithOffset(Matrix& M,
                                const Matrix& M2,
                                const size_t iOffset,
                                const size_t jOffset);
}

#endif //MATHUTIL_H
