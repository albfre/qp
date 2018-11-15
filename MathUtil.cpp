#include <iostream>
#include "MathUtil.h"

namespace MathUtil {
  extern "C" {
    // Solve linear system of equations
    // ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
    extern int dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);

    // LU factorize
    // (M, N, A, LDA, IPIV, INFO)
    extern int dgetrf_(int*, int*, double*, int*, int*, int*);

    // Solve Ax = b using LU factorization
    // (TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
    extern int dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
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
    int info = 1;

    // solve Ax = b
    dgesv_(&N, &nrhs, ArowMajor.data(), &N, p.data(), b.data(), &N, &info);

    if (info != 0) {
      std::cout << "Error: dgesv returned error code " << info << std::endl;
      return;
    }
  }

  void luFactorize(const Matrix& A, Vector& LU, std::vector<int>& piv) {
    const auto n = A.size();
    // Copy A to LU in row major format
    LU.resize(n * n);
    for (size_t i = 0; i < n; ++i) {
      assert(A[i].size() == n);
      for (size_t j = 0; j < n; ++j) {
        LU[j * n + i] = A[i][j];
      }
    }
    int N = static_cast<int>(n);
    int M = static_cast<int>(n);
    piv.resize(n);
    int info = 1;

    // factorize A
    dgetrf_(&M, &N, LU.data(), &M, piv.data(), &info);
    if (info != 0) {
      std::cout << "Error: dgetrf returned error code " << info << std::endl;
      return;
    }
  }

  void solveLinearSystemOfEquationsUsingLU(const Vector& LU, const std::vector<int>& piv, Vector& b) {
    const auto n = b.size();
    assert(piv.size() == n);
    assert(LU.size() == n * n);

    // allocate data
    char T = 'N'; // No transpose
    int N = static_cast<int>(n);
    int nrhs = 1;
    std::vector<int> p = piv;
    int info = 1;
    auto LUcopy = LU;

    // solve Ax = b
    dgetrs_(&T, &N, &nrhs, LUcopy.data(), &N, p.data(), b.data(), &N, &info);
    assert(LUcopy == LU);

    if (info != 0) {
      std::cout << "Error: dgetrs returned error code " << info << std::endl;
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
