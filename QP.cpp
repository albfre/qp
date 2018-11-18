#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>

#include "QP.h"
#include "LinearSolver.h"

namespace {
  using Vector = MathUtil::Vector;
  using Matrix = MathUtil::Matrix;

  // rA = Ax - y - b
  Vector computeRA(const Matrix& A, const Vector& x, const Vector& y, const Vector& b) {
    auto rA = MathUtil::matrixTimesVector(A, x);
    MathUtil::vectorMinusEqVector(rA, y);
    MathUtil::vectorMinusEqVector(rA, b);
    return rA;
  }

  // rQ = Qx - A' lambda + c
  Vector computeRQ(const Matrix& H, const Matrix& AT, const Vector& x, const Vector& lambda, const Vector& g) {
    auto rQ = MathUtil::matrixTimesVector(H, x);
    const auto ATLambda = MathUtil::matrixTimesVector(AT, lambda);
    MathUtil::vectorMinusEqVector(rQ, ATLambda);
    MathUtil::vectorPlusEqVector(rQ, g);
    return rQ;
  }

  // For matrix M = [H A'; A -Linv * Y], update -Linv * Y
  void updateMatrix(Matrix& M, const Vector& lambda, const Vector& y) {
    const auto size = M.size();
    const auto ySize = y.size();
    assert(lambda.size() == ySize);
    assert(size > ySize);
    const auto xSize = size - ySize;
    for (size_t i = 0; i < ySize; ++i) {
      const auto l = lambda[i];
      assert(l > 0.0);
      M[xSize + i][xSize + i] = -y[i] / l;
    }
  }

  auto computeSearchDirection(const Matrix& M,
                              const Vector& factorization,
                              const std::vector<int>& piv,
                              const Vector& rQ,
                              const Vector& rA,
                              const Matrix& A,
                              const Vector& y,
                              const Vector& lambda,
                              const double mu) {
    const auto xSize = rQ.size();
    Vector rhs(xSize);
    MathUtil::vectorMinusEqVector(rhs, rQ);
    Vector rhs2 = rA;
    MathUtil::vectorPlusEqVector(rhs2, y);
    MathUtil::vectorEqOpSelfVector(rhs2, lambda, [mu] (const auto& r, const auto& l) { assert(l > 0.0); return -r + mu / l; } );
    rhs.insert(rhs.end(), rhs2.cbegin(), rhs2.cend());

    // Compute dx, dlambda from M [dx; -dlambda] = [-rQ; -rA - y + mu./lambda]
    LinearSolver::solveLinearSystemOfEquationsSymmetricIndefinite(factorization, piv, rhs);
    Vector dx(rhs.begin(), rhs.begin() + xSize);
    Vector dlambda(rhs.begin() + xSize, rhs.end());
    std::transform(dlambda.cbegin(), dlambda.cend(), dlambda.begin(), std::negate<double>());

    // Compute dy = A dx + rA
    auto dy = MathUtil::matrixTimesVector(A, dx);
    MathUtil::vectorPlusEqVector(dy, rA);

    return std::make_tuple(dx, dy, dlambda);
  }

  double lineSearch(const Vector& y, const Vector& dy, const Vector& lambda, const Vector& dlambda) {
    const auto size = y.size();
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

double QP::objectiveValue(const Matrix& Q, const Vector& c, const Vector& x)
{
  const auto Qx = MathUtil::matrixTimesVector(Q, x);
  return 0.5 * MathUtil::dot(Qx, x) + MathUtil::dot(c, x);
}

QP::Solution QP::solveQP(const Matrix& Q,
                         const Vector& c,
                         const Matrix& A,
                         const Vector& b)
{
  Solution solution;

  // Preconditions
  const auto xSize = c.size();
  const auto ySize = b.size();
  assert(Q.size() == xSize);
  for (const auto& q : Q) {
    assert(q.size() == xSize);
  }
  assert(A.size() == ySize);
  for (const auto& a : A) {
    assert(a.size() == xSize);
  }

  // Handle unconstrained case
  if (ySize == 0) {
    auto x = c;
    LinearSolver::solveLinearSystemOfEquations(Q, x);
    solution.objectiveValue = objectiveValue(Q, c, x);
    solution.x = x;
    return solution;
  }

  // Store A'
  const auto AT = MathUtil::matrixTranspose(A);

  // Set up augmented system M = [H A'; A 0] (M22 is updated when computing search direction)
  const auto size = xSize + ySize;
  Matrix M(size, Vector(size));
  MathUtil::matrixEqMatrixWithOffset(M, Q, 0, 0);
  MathUtil::matrixEqMatrixWithOffset(M, AT, 0, xSize);
  MathUtil::matrixEqMatrixWithOffset(M, A, xSize, 0);

  // Initialize x, y, lambda
  Vector x(xSize, 1.0);
  // y = Ax - b
  auto y = computeRA(A, x, Vector(ySize), b);
  MathUtil::vectorEqOpSelfVector(y, y, [] (const auto& a, const auto& b) { return std::max(1e-3, b); });
  Vector lambda(ySize, 1.0);
  Vector prevX(xSize);
  Vector factorization(size * size);
  std::vector<int> piv(size);
  const auto fractionToBoundary = 0.99;

  for (size_t k = 0; k < 50; ++k ) {
    prevX = x;
    updateMatrix(M, lambda, y);
    LinearSolver::factorizeSymmetricIndefinite(M, factorization, piv);
    const auto rQ = computeRQ(Q, AT, x, lambda, c);
    const auto rA = computeRA(A, x, y, b);

    // Compute affine scaling step
    const auto [_, dyAff, dlambdaAff] = computeSearchDirection(M, factorization, piv, rQ, rA, A, y, lambda, 0.0);
    const auto alphaAff = lineSearch(y, dyAff, lambda, dlambdaAff);

    // Compute yAff, lambdaAff
    auto yAff = y;
    auto lambdaAff = lambda;
    MathUtil::vectorPlusEqScalarTimesVector(yAff, alphaAff, dyAff);
    MathUtil::vectorPlusEqScalarTimesVector(lambdaAff, alphaAff, dlambdaAff);
    const auto muAff = MathUtil::dot(yAff, lambdaAff) / ySize;

    // Compute centralizing step
    const auto mu = MathUtil::dot(y, lambda) / ySize;
    const auto sigma = std::pow(muAff / mu, 3.0); // Mehrotra heuristic
    const auto [dx, dy, dlambda] = computeSearchDirection(M, factorization, piv, rQ, rA, A, y, lambda, sigma * mu);
    const auto alpha = fractionToBoundary * lineSearch(y, dy, lambda, dlambda);

    // Update x, y, lambda
    MathUtil::vectorPlusEqScalarTimesVector(x, alpha, dx);
    MathUtil::vectorPlusEqScalarTimesVector(y, alpha, dy);
    MathUtil::vectorPlusEqScalarTimesVector(lambda, alpha, dlambda);

    std::cout << k << ". alphaAff: " << alphaAff << ", alpha: " << alpha << ", obj: " << objectiveValue(Q, c, x) << std::endl;
    const auto diff = std::inner_product(x.cbegin(), x.cend(), prevX.cbegin(), 0.0, std::plus<double>(),
                                         [](const auto& a, const auto& b) { return std::fabs(a - b); });
    if (diff < 1e-12) {
      break;
    }
  }

  solution.x = x;
  solution.objectiveValue = objectiveValue(Q, c, x);
  y = computeRA(A, x, Vector(ySize), b);
  solution.residual = std::accumulate(y.cbegin(), y.cend(), 0.0,
                                      [] (const auto& s, const auto& v) { return s + v < 0.0 ? std::fabs(v) : 0.0; });
  return solution;
}
