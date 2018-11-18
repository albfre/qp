#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <vector>

namespace LinearSolver {
  using Vector = std::vector<double>;
  using Matrix = std::vector<Vector>;

  void solveLinearSystemOfEquations(const Matrix& A, Vector& b);

  void factorizeGeneral(const Matrix& A, Vector& factorization, std::vector<int>& piv);
  void solveLinearSystemOfEquationsGeneral(const Vector& factorization, const std::vector<int>& piv, Vector& b);

  void factorizeSymmetricPositiveDefinite(const Matrix& A, Vector& factorization);
  void solveLinearSystemOfEquationsSymmetricPositiveDefinite(const Vector& factorization, Vector& b);

  void factorizeSymmetricIndefinite(const Matrix& A, Vector& factorization, std::vector<int>& piv);
  void solveLinearSystemOfEquationsSymmetricIndefinite(const Vector& factorization, const std::vector<int>& piv, Vector& b);

  void factorizeIndefinite(const Matrix& A, Vector& factorization, std::vector<int>& piv);
  void solveLinearSystemOfEquationsIndefinite(const Vector& factorization, const std::vector<int>& piv, Vector& b);
}


#endif //LINEAR_SOLVER_H
