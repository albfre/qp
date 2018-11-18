#ifndef QP_H
#define QP_H

#include "MathUtil.h"

namespace QP {
  using Vector = MathUtil::Vector;
  using Matrix = MathUtil::Matrix;

  struct Solution {
    double objectiveValue = 1e50;
    double residual = 0.0;
    Vector x;
  };

  double objectiveValue(const Matrix& Q,
                        const Vector& c,
                        const Vector& x);

  /* Solve problem
   * min 0.5 x' Q x + c' x
   * st  Ax >= b
   */
  Solution solveQP(const Matrix& Q,
                   const Vector& c,
                   const Matrix& A,
                   const Vector& b);
}

#endif //QP_H
