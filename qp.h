#ifndef QP_H
#define QP_H

#include <algorithm>
#include <memory>
#include <vector>

namespace QP {
  struct Solution {
    double objectiveValue;
    std::vector<double> x;
  };

  using Matrix = std::vector<std::vector<double>>;
  using Vector = std::vector<double>;

  double objectiveValue(const Matrix& H,
                        const Vector& g,
                        const Vector& x);

  /* Solve problem
   * min 0.5 x' H x + g' x
   * st  Ax >= b
   */
  Solution solveQP(const Matrix& H,
                   const Vector& g,
                   const Matrix& A,
                   const Vector& b);
}


#endif //QP_H
