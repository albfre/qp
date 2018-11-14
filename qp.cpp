#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>
#include <cmath>
#include <numeric>
#include <functional>
#include <sstream>
#include <iostream>
#include <tuple>

#include "qp.h"

namespace {
  extern "C" {
    //              ( N, NRHS, A,       LDA,     IPIV,    B,     LDB,  INFO)
    extern int dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
  }

  using Matrix = QP::Matrix;
  using Vector = QP::Vector;

  void print(const Vector& x) {
    std::stringstream ss;
    for (const auto& xi : x) {
      ss << xi << ", ";
    }
    ss << std::endl;
    std::cout << ss.str();
  }

  void solveLinearSystemOfEquations_(const std::vector<std::vector<double>>& A, std::vector<double>& b)
  {
    const auto n = b.size();
    assert(A.size() == n);
    std::vector<double> ArowMajor(n * n);
    for (size_t i = 0; i < n; ++i) {
      assert(A[i].size() == n);
      for (size_t j = 0; j < n; ++j) {
        //std::cout << A[i][j] << " ";
        ArowMajor[j * n + i] = A[i][j];
      }
      //std::cout << std::endl;
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

  Matrix matrixTranspose(const Matrix& m, const size_t xSize, const size_t ySize)
  {
    assert(m.size() == ySize);
    for (const auto& row : m) {
      assert(row.size() == xSize);
    }
    Matrix mT(xSize, Vector(ySize));
    for (size_t i = 0; i < ySize; ++i) {
      for (size_t j = 0; j < xSize; ++j) {
        mT[j][i] = m[i][j];
      }
    }
    return mT;
  }

  template<typename Op>
  void vectorEqOpSelfVector(Vector& r, const Vector& v, const Op& op) {
    const auto size = r.size();
    assert(v.size() == size);
    for (size_t i = 0; i < size; ++i) {
      r[i] = op(r[i], v[i]);
    }
  }

  template<class Op>
  void matrixEqOpMatrixWithOffset(Matrix& M, const Matrix& M2, const size_t iOffset, const size_t jOffset, const Op& op) {
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

  // rp = Ax - y - b
  Vector computeRp(const Matrix& A, const Vector& x, const Vector& y, const Vector& b) {
    const auto ySize = y.size();
    Vector rp(ySize);
    matrixTimesVector(A, x, rp);
    vectorEqOpSelfVector(rp, y, std::minus<double>());
    vectorEqOpSelfVector(rp, b, std::minus<double>());
    return rp;
  }

  // rd = Hx - A' lambda + g
  Vector computeRd(const Matrix& H, const Matrix& AT, const Vector& x, const Vector& lambda, const Vector& g) {
    const auto xSize = x.size();
    Vector rd(xSize);
    matrixTimesVector(H, x, rd);
    Vector ATLambda(xSize);
    matrixTimesVector(AT, lambda, ATLambda);
    vectorEqOpSelfVector(rd, ATLambda, std::minus<double>());
    vectorEqOpSelfVector(rd, g, std::plus<double>());
    return rd;
  }

  // For matrix M = [H A'; A Linv * Y], update Linv * Y
  void updateMatrix(Matrix& M, const Vector& lambda, const Vector& y) {
    const auto size = M.size();
    const auto ySize = y.size();
    assert(lambda.size() == ySize);
    assert(size > ySize);
    const auto xSize = size - ySize;
    for (size_t i = 0; i < ySize; ++i) {
      const auto l = lambda[i];
      assert(l > 0.0);
      M[xSize + i][xSize + i] = y[i] / l;
    }
  }

  auto computeSearchDirection(Matrix& M,
                              const Matrix& H,
                              const Vector& g,
                              const Matrix& A,
                              const Matrix& AT,
                              const Vector& b,
                              const Vector& x,
                              const Vector& y,
                              const Vector& lambda,
                              const double mu) {
    const auto xSize = x.size();
    const auto ySize = y.size();
    const auto size = xSize + ySize;

    updateMatrix(M, lambda, y);
    const auto rd = computeRd(H, AT, x, lambda, g);
    const auto rp = computeRp(A, x, y, b);
    auto rpPlusY = rp;
    vectorEqOpSelfVector(rpPlusY, y, std::plus<double>());

    Vector rhs(size);
    for (size_t i = 0; i < xSize; ++i) {
      rhs[i] = -rd[i];
    }
    for (size_t i = 0; i < ySize; ++i) {
      rhs[xSize + i] = -rpPlusY[i] + mu / lambda[i];
    }

    // Compute dx, dlambda from M [dx; dlambda] = [-rd; -rp - y + mu./lambda]
    solveLinearSystemOfEquations_(M, rhs);
    Vector dx(rhs.begin(), rhs.begin() + xSize);
    Vector dlambda(rhs.begin() + xSize, rhs.end());

    // Compute dy = A dx + rp
    Vector dy(ySize);
    matrixTimesVector(A, dx, dy);
    vectorEqOpSelfVector(dy, rp, std::plus<double>());

    return std::make_tuple(dx, dy, dlambda);
  }

  double lineSearch(const Vector& y, const Vector& dy, const Vector& lambda, const Vector& dlambda) {
    auto size = y.size();
    assert(dy.size() == size);
    assert(lambda.size() == size);
    assert(dlambda.size() == size);

    auto alpha = 1.0;
    for (size_t i = 0; i < y.size(); ++i) {
      if (dy[i] < 0.0) {
        alpha = std::min(alpha, -y[i]/dy[i]);
      }
      if (dlambda[i] < 0.0) {
        alpha = std::min(alpha, -lambda[i]/dlambda[i]);
      }
    }
    if (alpha < 1e-15) {
      assert( false );
    }
    return alpha;
  }
}

double QP::objectiveValue(const Matrix& H, const Vector& g, const Vector& x)
{
  Vector Hx(x.size());
  matrixTimesVector(H, x, Hx);
  return 0.5 * dot(Hx, x) + dot(g, x);
}

QP::Solution QP::solveQP(const Matrix& H,
                         const Vector& g,
                         const Matrix& A,
                         const Vector& b)
{
  Solution solution;

  // Preconditions
  const auto xSize = g.size();
  const auto ySize = b.size();
  assert(H.size() == xSize);
  for (const auto& h : H) {
    assert(h.size() == xSize);
  }
  for (const auto& a: A) {
    assert(a.size() == xSize);
  }
  assert(A.size() == ySize);

  // Handle unconstrained case
  if (ySize == 0) {
    auto x = g;
    solveLinearSystemOfEquations_(H, x);
    solution.objectiveValue = objectiveValue(H, g, x);
    solution.x = x;
    return solution;
  }

  // Store A'
  const auto AT = matrixTranspose(A, xSize, ySize);

  // Set up augmented system M = [H A'; A 0] (M22 is updated when computing search direction)
  const auto size = xSize + ySize;
  Matrix M(size, Vector(size));
  const auto identity = [](const auto& a) { return a; };
  const auto negate = [](const auto& a) { return -a; };
  matrixEqOpMatrixWithOffset(M, H, 0, 0, identity);
  matrixEqOpMatrixWithOffset(M, AT, 0, xSize, negate);
  matrixEqOpMatrixWithOffset(M, A, xSize, 0, identity);

  // Initialize x, y, lambda
  Vector x(xSize, 1.0);
  // y = Ax - b
  Vector Ax(ySize);
  matrixTimesVector(A, x, Ax);
  auto y = Ax;
  vectorEqOpSelfVector(y, b, std::minus<double>());
  vectorEqOpSelfVector(y, y, [] (const auto& a, const auto& b) { return std::max(1e-3, b); });
  Vector lambda(ySize, 1.0);
  Vector yAff(ySize);
  Vector lambdaAff(ySize);
  Vector prevX(xSize);

  for (size_t k = 0; k < 50; ++k ) {
    prevX = x;

    // Compute affine scaling step
    yAff = y;
    lambdaAff = lambda;
    const auto [_, dyAff, dlambdaAff] = computeSearchDirection(M, H, g, A, AT, b, x, y, lambda, 0.0);
    const auto alphaAff = lineSearch(y, dyAff, lambda, dlambdaAff);

    // Update yAff, lambdaAff
    const auto scalarAddAff = [alphaAff] (const auto& r, const auto& v) { return r + alphaAff * v; };
    vectorEqOpSelfVector(yAff, dyAff, scalarAddAff);
    vectorEqOpSelfVector(lambdaAff, dlambdaAff, scalarAddAff );
    const auto muAff = dot(yAff, lambdaAff) / ySize;

    // Compute centralizing step
    const auto mu = dot(y, lambda) / ySize;
    const auto div = muAff / mu;
    const auto sigma = div * div * div;
    const auto [dx, dy, dlambda] = computeSearchDirection(M, H, g, A, AT, b, x, y, lambda, sigma * mu);
    const auto alpha = lineSearch(y, dy, lambda, dlambda);

    // Update x, y, lambda
    const auto scalarAdd = [alpha] (const auto& r, const auto& v) { return r + alpha * v; };
    vectorEqOpSelfVector(x, dx, scalarAdd);
    vectorEqOpSelfVector(y, dy, scalarAdd);
    vectorEqOpSelfVector(lambda, dlambda, scalarAdd);

    std::cout << k << ". alphaAff: " << alphaAff << ", alpha: " << alpha << ", obj: " << objectiveValue(H, g, x) << std::endl;
    const auto diff = std::inner_product(x.cbegin(), x.cend(), prevX.cbegin(), 0.0, std::plus<double>(),
                                         [](const auto& a, const auto& b) { return std::fabs(a - b); });
    if (diff < 1e-12) {
      break;
    }
  }

  solution.x = x;
  solution.objectiveValue = objectiveValue(H, g, x);
  return solution;
}
