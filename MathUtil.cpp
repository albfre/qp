#include <iostream>
#include "MathUtil.h"

namespace MathUtil {
  extern "C" {
    //              ( N, NRHS, A,       LDA,     IPIV,    B,     LDB,  INFO)
    extern int dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
  }

  void solveLinearSystemOfEquations(const Matrix& A, Vector& b)
  {
    const auto n = b.size();
    assert(A.size() == n);
    std::vector<double> ArowMajor(n * n);
    for (size_t i = 0; i < n; ++i) {
      assert(A[i].size() == n);
      for (size_t j = 0; j < n; ++j) {
        ArowMajor[j * n + i] = A[i][j];
      }
    }
    // allocate data
    int N = static_cast<int>(n);
    int nrhs = 1;
    std::vector<int> p(n);
    int info;

    // solve Ax = b
    dgesv_(&N, &nrhs, ArowMajor.data(), &N, p.data(), b.data(), &N, &info);

    if (info != 0) {
      std::cout << "Error: dgesv returned error code " << info << std::endl;
      return;
    }
  }

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

  void vectorPlusEqVector(Vector& r, const Vector& v)
  {
    vectorEqOpSelfVector(r, v, std::plus<double>());
  }

  void vectorMinusEqVector(Vector& r, const Vector& v) 
  {
    vectorEqOpSelfVector(r, v, std::minus<double>());
  }

  void vectorPlusEqScalarTimesVector(Vector& r, const double scalar, const Vector& v)
  {
    vectorEqOpSelfVector(r, v, [scalar] (const auto& rx, const auto& vx) { return rx + scalar * vx; });
  }
}
