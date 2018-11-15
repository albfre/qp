#include <iostream>
#include <cmath>
#include <tuple>

#include "QP.h"

namespace {
  using Vector = MathUtil::Vector;
  using Matrix = MathUtil::Matrix;

  // rp = Ax - y - b
  Vector computeRp(const Matrix& A, const Vector& x, const Vector& y, const Vector& b) {
    auto rp = MathUtil::matrixTimesVector(A, x);
    MathUtil::vectorMinusEqVector(rp, y);
    MathUtil::vectorMinusEqVector(rp, b);
    return rp;
  }

  // rd = Hx - A' lambda + g
  Vector computeRd(const Matrix& H, const Matrix& AT, const Vector& x, const Vector& lambda, const Vector& g) {
    auto rd = MathUtil::matrixTimesVector(H, x);
    const auto ATLambda = MathUtil::matrixTimesVector(AT, lambda);
    MathUtil::vectorMinusEqVector(rd, ATLambda);
    MathUtil::vectorPlusEqVector(rd, g);
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

  auto computeSearchDirection(const Matrix& M,
                              const Vector& LU,
                              const std::vector<int>& piv,
                              const Vector& rd,
                              const Vector& rp,
                              const Matrix& A,
                              const Vector& y,
                              const Vector& lambda,
                              const double mu) {
    const auto xSize = rd.size();
    Vector rhs(xSize);
    MathUtil::vectorMinusEqVector(rhs, rd);
    Vector rhs2 = rp;
    MathUtil::vectorPlusEqVector(rhs2, y);
    MathUtil::vectorEqOpSelfVector(rhs2, lambda, [mu] (const auto& r, const auto& l) { assert(l > 0.0); return -r + mu / l; } );
    rhs.insert(rhs.end(), rhs2.cbegin(), rhs2.cend());

    // Compute dx, dlambda from M [dx; dlambda] = [-rd; -rp - y + mu./lambda]
    MathUtil::solveLinearSystemOfEquationsUsingLU(LU, piv, rhs);
    Vector dx(rhs.begin(), rhs.begin() + xSize);
    Vector dlambda(rhs.begin() + xSize, rhs.end());

    // Compute dy = A dx + rp
    auto dy = MathUtil::matrixTimesVector(A, dx);
    MathUtil::vectorPlusEqVector(dy, rp);

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

double QP::objectiveValue(const Matrix& H, const Vector& g, const Vector& x)
{
  const auto Hx = MathUtil::matrixTimesVector(H, x);
  return 0.5 * MathUtil::dot(Hx, x) + MathUtil::dot(g, x);
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
    MathUtil::solveLinearSystemOfEquations(H, x);
    solution.objectiveValue = objectiveValue(H, g, x);
    solution.x = x;
    return solution;
  }

  // Store A'
  const auto AT = MathUtil::matrixTranspose(A);

  // Set up augmented system M = [H A'; A 0] (M22 is updated when computing search direction)
  const auto size = xSize + ySize;
  Matrix M(size, Vector(size));
  const auto identity = [](const auto& a) { return a; };
  const auto negate = [](const auto& a) { return -a; };
  MathUtil::matrixEqOpMatrixWithOffset(M, H, 0, 0, identity);
  MathUtil::matrixEqOpMatrixWithOffset(M, AT, 0, xSize, negate);
  MathUtil::matrixEqOpMatrixWithOffset(M, A, xSize, 0, identity);

  // Initialize x, y, lambda
  Vector x(xSize, 1.0);
  // y = Ax - b
  auto y = computeRp(A, x, Vector(ySize), b);
  MathUtil::vectorEqOpSelfVector(y, y, [] (const auto& a, const auto& b) { return std::max(1e-3, b); });
  Vector lambda(ySize, 1.0);
  Vector prevX(xSize);
  Vector LU(size * size);
  std::vector<int> piv(size);
  const auto fractionToBoundary = 0.99;

  for (size_t k = 0; k < 50; ++k ) {
    prevX = x;
    updateMatrix(M, lambda, y);
    MathUtil::luFactorize(M, LU, piv);
    const auto rd = computeRd(H, AT, x, lambda, g);
    const auto rp = computeRp(A, x, y, b);

    // Compute affine scaling step
    const auto [_, dyAff, dlambdaAff] = computeSearchDirection(M, LU, piv, rd, rp, A, y, lambda, 0.0);
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
    const auto [dx, dy, dlambda] = computeSearchDirection(M, LU, piv, rd, rp, A, y, lambda, sigma * mu);
    const auto alpha = fractionToBoundary * lineSearch(y, dy, lambda, dlambda);

    // Update x, y, lambda
    MathUtil::vectorPlusEqScalarTimesVector(x, alpha, dx);
    MathUtil::vectorPlusEqScalarTimesVector(y, alpha, dy);
    MathUtil::vectorPlusEqScalarTimesVector(lambda, alpha, dlambda);

    std::cout << k << ". alphaAff: " << alphaAff << ", alpha: " << alpha << ", obj: " << objectiveValue(H, g, x) << std::endl;
    const auto diff = std::inner_product(x.cbegin(), x.cend(), prevX.cbegin(), 0.0, std::plus<double>(),
                                         [](const auto& a, const auto& b) { return std::fabs(a - b); });
    if (diff < 1e-13) {
      break;
    }
  }

  solution.x = x;
  solution.objectiveValue = objectiveValue(H, g, x);
  return solution;
}
