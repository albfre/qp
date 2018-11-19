#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>

#include "QP.h"
#include "LinearSolver.h"

namespace {
  using Vector = MathUtil::Vector;
  using Matrix = MathUtil::Matrix;

  // rQ = Qx + c - A' y - C' z
  Vector computeRQ(const Matrix& Q, const Matrix& AT, const Matrix& CT, const Vector& c, const Vector& x, const Vector& y, const Vector& z) {
    auto rQ = MathUtil::matrixTimesVector(Q, x);
    MathUtil::vectorPlusEqVector(rQ, c);
    if (!y.empty()) {
      const auto ATy = MathUtil::matrixTimesVector(AT, y);
      MathUtil::vectorMinusEqVector(rQ, ATy);
    }
    if (!z.empty()) {
      const auto CTz = MathUtil::matrixTimesVector(CT, z);
      MathUtil::vectorMinusEqVector(rQ, CTz);
    }
    return rQ;
  }

  // rA = Ax - b
  Vector computeRA(const Matrix& A, const Vector& b, const Vector& x) {
    auto rA = MathUtil::matrixTimesVector(A, x);
    MathUtil::vectorMinusEqVector(rA, b);
    return rA;
  }

  // rC = Cx - s - d
  Vector computeRC(const Matrix& C, const Vector& d, const Vector& x, const Vector& s) {
    auto rC = MathUtil::matrixTimesVector(C, x);
    MathUtil::vectorMinusEqVector(rC, s);
    MathUtil::vectorMinusEqVector(rC, d);
    return rC;
  }

  // rzs = ZSe - sigma * mu e
  Vector computeRzs(const Vector& z, const Vector& s, const double sigmaTimesMu) {
    auto rzs = z;
    MathUtil::vectorEqOpSelfVector(rzs, s, [sigmaTimesMu] (const auto& zi, const auto& si) { return zi * si - sigmaTimesMu; });
    return rzs;
  }

  double computeInfeasibility(const Matrix& A, const Matrix& C, const Vector& b, const Vector& d, const Vector& x) {
    const auto rA = computeRA(A, b, x);
    auto r = std::accumulate(rA.cbegin(), rA.cend(), 0.0,
                             [] (const auto& s, const auto& v) { return s + std::fabs(v); });
    const auto rC = computeRC(C, d, x, Vector(C.size()));
    r += std::accumulate(rC.cbegin(), rC.cend(), 0.0,
                         [] (const auto& s, const auto& v) { return s + (v < 0.0 ? std::fabs(v) : 0.0); });
    return r;
  }

  // For matrix M = [Q A' C'; A 0 0; C 0 -Zinv * S], update -Zinv * S
  void updateMatrix(Matrix& M, const Vector& s, const Vector& z) {
    const auto size = M.size();
    const auto nz = z.size();
    assert(s.size() == nz);
    assert(size > nz);
    const auto nxPlusny = size - nz;
    for (size_t i = 0; i < nz; ++i) {
      const auto zi = z[i];
      assert(zi > 0.0);
      M[nxPlusny + i][nxPlusny + i] = -s[i] / zi;
    }
  }

 /** Compute search direction by solving the system
   *
   *  [ Q A' C'      ] [  dx ]    [ rQ = Qx + c + A'y - C'z               ]
   *  [ A 0  0       ] [ -dy ] = -[ rA = Ax - b                           ]
   *  [ C 0  -Zinv S ] [ -dz ]    [ rC + Zinv rzs = Cx - s - d + Zinv rzs ]
   *
   *  ds = -Zinv (rzs + S dz)
   *
   *  where rzs = ZSe - sigma mu e
   *  i.e., rhs3 = -rC - s + sigma * mu ./ z
   */
  auto computeSearchDirection(const Matrix& M,
                              const Matrix& C,
                              const Vector& factorization,
                              const std::vector<int>& piv,
                              const Vector& rQ,
                              const Vector& rA,
                              const Vector& rC,
                              const Vector& rzs,
                              const Vector& y,
                              const Vector& z,
                              const Vector& s) {
    const auto nx = rQ.size();
    const auto ny = rA.size();
    Vector rhs(nx);
    MathUtil::vectorMinusEqVector(rhs, rQ); // rhs1 = -rQ

    Vector rhs2(ny);
    MathUtil::vectorMinusEqVector(rhs2, rA); // rhs2 = -rA
    rhs.insert(rhs.end(), rhs2.cbegin(), rhs2.cend());

    Vector rhs3 = rzs;
    MathUtil::vectorDivEqMinusVector(rhs3, z); // rhs3 = -Zinv rzs
    MathUtil::vectorMinusEqVector(rhs3, rC); // rhs3 = -rC - Zinv rzs
    rhs.insert(rhs.end(), rhs3.cbegin(), rhs3.cend());

    // Compute dx, dy, dz from M [dx; -dy; -dz] = [-rQ; -rA; -rC - Zinv rzs]
    LinearSolver::solveLinearSystemOfEquationsSymmetricIndefinite(factorization, piv, rhs);
    Vector dx(rhs.begin(), rhs.begin() + nx);
    Vector dy(rhs.begin() + nx, rhs.begin() + nx + ny);
    Vector dz(rhs.begin() + nx + ny, rhs.end());
    std::transform(dy.cbegin(), dy.cend(), dy.begin(), std::negate<double>());
    std::transform(dz.cbegin(), dz.cend(), dz.begin(), std::negate<double>());

    // Compute ds = -Zinv (rzs + S dz)
    // Compute ds = C dx + rC
//    auto ds = dz;
//    MathUtil::vectorTimesEqVector(ds, s);
//    MathUtil::vectorPlusEqVector(dz, rzs);
//    MathUtil::vectorDivEqMinusVector(dz, z);
    auto ds = MathUtil::matrixTimesVector(C, dx);
    MathUtil::vectorPlusEqVector(ds, rC);

    return std::make_tuple(dx, dy, dz, ds);
  }

  double lineSearch(const Vector& z, const Vector& dz, const Vector& s, const Vector& ds) {
    const auto size = z.size();
    assert(dz.size() == size);
    assert(s.size() == size);
    assert(ds.size() == size);

    auto alpha = 1.0;
    for (size_t i = 0; i < size; ++i) {
      if (dz[i] < 0.0) {
        alpha = std::min(alpha, -z[i]/dz[i]);
      }
      if (ds[i] < 0.0) {
        alpha = std::min(alpha, -s[i]/ds[i]);
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
                         const Vector& b,
                         const Matrix& C,
                         const Vector& d)
{
  // Preconditions
  const auto nx = c.size();
  const auto ny = b.size();
  const auto nz = d.size(); 
  const auto assertMatrix = [] (const auto& M, const auto& m, const auto& n) {
    assert(M.size() == m);
    assert(std::all_of(M.cbegin(), M.cend(), [n] (const auto& m) { return m.size() == n; }));
  };
  assertMatrix(Q, nx, nx);
  assertMatrix(A, ny, nx);
  assertMatrix(C, nz, nx);

  // Store A' and C'
  const auto AT = MathUtil::matrixTranspose(A);
  const auto CT = MathUtil::matrixTranspose(C);

  // Set up augmented system M = [Q A' C'; A 0 0; C 0 0] (M33 is updated when computing search direction)
  const auto size = nx + ny + nz;
  Matrix M(size, Vector(size));
  MathUtil::matrixEqMatrixWithOffset(M, Q, 0, 0);
  MathUtil::matrixEqMatrixWithOffset(M, AT, 0, nx);
  MathUtil::matrixEqMatrixWithOffset(M, CT, 0, nx + ny);
  MathUtil::matrixEqMatrixWithOffset(M, A, nx, 0);
  MathUtil::matrixEqMatrixWithOffset(M, C, nx + ny, 0);

  // Handle unconstrained case
  if (nz == 0) {
    Vector rhs(c.size());
    MathUtil::vectorMinusEqVector(rhs, c);
    rhs.insert(rhs.end(), b.cbegin(), b.cend());
    LinearSolver::solveLinearSystemOfEquations(M, rhs);
    const Vector x(rhs.cbegin(), rhs.cbegin() + nx);
    return Solution{ x, objectiveValue(Q, c, x), computeInfeasibility(A, C, b, d, x) };
  }

  // Initialize x, y, z, s
  Vector x(nx, 1.0);
  Vector y(ny, 1.0);
  Vector z(nz, 1.0);
  auto s = computeRC(C, d, x, Vector(nz));
  MathUtil::vectorEqOpSelfVector(s, s, [] (const auto& a, const auto& d) { return std::max(1e-3, d); });

  Vector factorization;
  std::vector<int> piv(size);
  const auto fractionToBoundary = 0.99;

  for (size_t k = 0; k < 100; ++k ) {
    updateMatrix(M, s, z);
    LinearSolver::factorizeSymmetricIndefinite(M, factorization, piv);
    const auto rQ = computeRQ(Q, AT, CT, c, x, y, z);
    const auto rA = computeRA(A, b, x);
    const auto rC = computeRC(C, d, x, s);

    // Compute affine scaling step
    const auto rzsAff = computeRzs(z, s, 0.0);
    const auto [dxAff, dyAff, dzAff, dsAff] = computeSearchDirection(M, C, factorization, piv, rQ, rA, rC, rzsAff, y, z, s);
    const auto alphaAff = lineSearch(z, dzAff, s, dsAff);

    // Compute zAff, sAff
    auto zAff = z;
    auto sAff = s;
    MathUtil::vectorPlusEqScalarTimesVector(zAff, alphaAff, dzAff);
    MathUtil::vectorPlusEqScalarTimesVector(sAff, alphaAff, dsAff);
    const auto muAff = MathUtil::dot(zAff, sAff) / nz;

    // Compute aggregated centering-corrector direction
    const auto mu = MathUtil::dot(z, s) / nz;
    const auto sigma = std::pow(muAff / mu, 3.0);
    const auto rzsCorr = computeRzs(dzAff, dsAff, 0.0);
    auto rzsCentCorr = computeRzs(z, s, sigma * mu);
    MathUtil::vectorPlusEqVector(rzsCentCorr, rzsCorr);


    const auto [dx, dy, dz, ds] = computeSearchDirection(M, C, factorization, piv, rQ, rA, rC, rzsCentCorr, y, z, s);
    const auto alpha = fractionToBoundary * lineSearch(z, dz, s, ds);

    // Update x, y, z, s
    MathUtil::vectorPlusEqScalarTimesVector(x, alpha, dx);
    MathUtil::vectorPlusEqScalarTimesVector(y, alpha, dy);
    MathUtil::vectorPlusEqScalarTimesVector(z, alpha, dz);
    MathUtil::vectorPlusEqScalarTimesVector(s, alpha, ds);
    const auto gap = MathUtil::dot(z, s) / nz;

    std::cout << k << ". alphaAff: " << alphaAff << ", alpha: " << alpha << ", gap: " << gap << ", obj: " << objectiveValue(Q, c, x) << std::endl;
    if (gap < 1e-12) {
      break;
    }
  }

  return Solution{ x, objectiveValue(Q, c, x), computeInfeasibility(A, C, b, d, x) };
}
