#include <iostream>
#include <cmath>
#include <cassert>
#include <numeric>
#include "LinearSolver.h"

namespace {
  using Vector = LinearSolver::Vector;
  using Matrix = LinearSolver::Matrix;

  void matrixToRowMajorFormat_(const Matrix& A, Vector& ArowMajor) {
    const auto n = A.size();
    // Copy A to LU in row major format
    ArowMajor.resize(n * n);
    for (size_t i = 0; i < n; ++i) {
      assert(A[i].size() == n);
      for (size_t j = 0; j < n; ++j) {
        ArowMajor[j * n + i] = A[i][j];
      }
    }
  }

  void matrixToPackedLowerFormat_(const Matrix& A, Vector& Apacked) {
    const auto n = A.size();
    // Copy A to factorization in packed lower Lapack format
    Apacked.resize(n * (n + 1) / 2);
    size_t fi = 0;
    for (size_t i = 0; i < n; ++i) {
      assert(A[i].size() == n);
      for (size_t j = i; j < n; ++j, ++fi) {
        Apacked[fi] = A[i][j];
      }
    }
    assert(fi == Apacked.size());
  }

  void handleErrors_(const int info, const std::string& functionName) {
    if (info != 0) {
      std::cout << "Error: " + functionName + " returned error code " << info << std::endl;
    }
  }
}

namespace LinearSolver {
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

    // Cholesky factorize
    // (UPLO, N, A, LDA, INFO)
    extern int dpotrf_(char*, int*, double*, int*, int*);

    // Solve Ax = b using Cholesky factorization
    // (UPLO, N, NRHS, A, LDA, B, LDB, INFO)
    extern int dpotrs_(char*, int*, int*, double*, int*, double*, int*, int*);

    // Factorize symmetric indefinite matrix using Bunch-Kaufman diagonal pivoting
    // (UPLO, N, AP, IPIV, INFO)
    extern int dsptrf_(char*, int*, double*, int*,  int*);

    // Solve Ax = b using symmetric indefinite factorization, A packed
    // (UPLO, N, NRHS, AP, IPIV, B, LDB, INFO)
    extern int dsptrs_(char*, int*, int*, double*, int*, double*, int*, int*);

    // Factorize indefinite matrix
    // (UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO)
    extern int dsytrf_(char*, int*, double*, int*, int*, double*, int*,  int*);

    // Solve Ax = b using indefinite factorization
    // (UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
    extern int dsytrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
  }

  void solveLinearSystemOfEquations(const Matrix& A, Vector& b)
  {
    Vector ArowMajor;
    matrixToRowMajorFormat_(A, ArowMajor);
    const auto n = A.size();
    std::vector<int> p(n);

    auto N = static_cast<int>(n);
    auto nrhs = 1;
    auto info = 1;

    // solve Ax = b by LU factorization
    dgesv_(&N, &nrhs, ArowMajor.data(), &N, p.data(), b.data(), &N, &info);
    handleErrors_(info, "dgesv");
  }

  void factorizeGeneral(const Matrix& A, Vector& factorization, std::vector<int>& piv) {
    matrixToRowMajorFormat_(A, factorization);
    const auto n = A.size();
    piv.resize(n);

    auto N = static_cast<int>(n);
    auto M = static_cast<int>(n);
    auto info = 1;

    // LU factorize A
    dgetrf_(&M, &N, factorization.data(), &M, piv.data(), &info);
    handleErrors_(info, "dgetrf");
  }

  void solveLinearSystemOfEquationsGeneral(const Vector& factorization, const std::vector<int>& piv, Vector& b) {
    const auto n = b.size();
    assert(piv.size() == n);
    assert(factorization.size() == n * n);

    auto T = 'N'; // No transpose
    auto N = static_cast<int>(n);
    auto nrhs = 1;
    auto p = piv;
    auto info = 1;

    // solve Ax = b using LU factorization
    auto factorizationCopy = factorization;
    dgetrs_(&T, &N, &nrhs, factorizationCopy.data(), &N, p.data(), b.data(), &N, &info);
    assert(factorizationCopy == factorization);
    handleErrors_(info, "dgetrs");
  }

  void factorizeSymmetricPositiveDefinite(const Matrix& A, Vector& factorization) {
    matrixToRowMajorFormat_(A, factorization);
    const auto n = A.size();

    auto UPLO = 'U';
    auto N = static_cast<int>(n);
    auto info = 1;

    // Cholesky factorize A
    dpotrf_(&UPLO, &N, factorization.data(), &N, &info);
    handleErrors_(info, "dpotrf");
  }

  void solveLinearSystemOfEquationsSymmetricPositiveDefinite(const Vector& factorization, Vector& b) {
    const auto n = b.size();
    assert(factorization.size() == n * n);

    // allocate data
    auto UPLO = 'U';
    auto N = static_cast<int>(n);
    auto nrhs = 1;
    auto info = 1;

    // solve Ax = b using Cholesky factorization
    auto factorizationCopy = factorization;
    dpotrs_(&UPLO, &N, &nrhs, factorizationCopy.data(), &N, b.data(), &N, &info);
    assert(factorizationCopy == factorization);
    handleErrors_(info, "dpotrs");
  }

  void factorizeSymmetricIndefinite(const Matrix& A, Vector& factorization, std::vector<int>& piv) {
    matrixToPackedLowerFormat_(A, factorization);
    const auto n = A.size();
    piv.resize(n);

    auto UPLO = 'L';
    auto N = static_cast<int>(n);
    auto info = 1;

    // factorize A using Bunch-Kaufman diagonal pivoting
    dsptrf_(&UPLO, &N, factorization.data(), piv.data(), &info);
    handleErrors_(info, "dsptrf");
  }

  void solveLinearSystemOfEquationsSymmetricIndefinite(const Vector& factorization, const std::vector<int>& piv, Vector& b) {
    const auto n = b.size();
    assert(piv.size() == n);
    assert(factorization.size() == n * (n + 1) / 2);

    // allocate data
    auto UPLO = 'L';
    auto N = static_cast<int>(n);
    auto nrhs = 1;
    auto p = piv;
    auto info = 1;

    // solve Ax = b
    // (UPLO, N, NRHS, A, IPIV, B, LDB, INFO)
    auto factorizationCopy = factorization;
    dsptrs_(&UPLO, &N, &nrhs, factorizationCopy.data(), p.data(), b.data(), &N, &info);
    assert(factorizationCopy == factorization);
    handleErrors_(info, "dsptrs");
  }

  void factorizeIndefinite(const Matrix& A, Vector& factorization, std::vector<int>& piv) {
    matrixToRowMajorFormat_(A, factorization);
    const auto n = A.size();
    piv.resize(n);

    auto UPLO = 'L';
    auto N = static_cast<int>(n);
    auto info = 1;
    const auto nb = 10;
    std::vector<double> WORK(n * nb);
    auto LWORK = static_cast<int>(WORK.size());

    // factorize A
    dsytrf_(&UPLO, &N, factorization.data(), &N, piv.data(), WORK.data(), &LWORK, &info);
    handleErrors_(info, "dsytrf");
  }

  void solveLinearSystemOfEquationsIndefinite(const Vector& factorization, const std::vector<int>& piv, Vector& b) {
    const auto n = b.size();
    assert(piv.size() == n);
    assert(factorization.size() == n * n);

    auto UPLO = 'L';
    auto N = static_cast<int>(n);
    auto nrhs = 1;
    auto p = piv;
    auto info = 1;

    // solve Ax = b
    auto factorizationCopy = factorization;
    dsytrs_(&UPLO, &N, &nrhs, factorizationCopy.data(), &N, p.data(), b.data(), &N, &info);
    assert(factorizationCopy == factorization);
    handleErrors_(info, "dsytrs");
  }
}
