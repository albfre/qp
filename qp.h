#ifndef QP_H
#define QP_H

#include "MathUtil.h"

namespace QP {
  using Vector = MathUtil::Vector;
  using Matrix = MathUtil::Matrix;

  struct Solution {
    double objectiveValue;
    Vector x;
  };

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
