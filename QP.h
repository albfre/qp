#ifndef QP_H
#define QP_H

#include "MathUtil.h"

namespace QP {
  using Vector = MathUtil::Vector;
  using Matrix = MathUtil::Matrix;

  struct Solution {
    Vector x;
    double objectiveValue;
    double infeasibility;
  };

  double objectiveValue(const Matrix& Q,
                        const Vector& c,
                        const Vector& x);

  /**
   *  Solve problem
   *  min 0.5 x'Qx + c'x
   *  st  Ax == b
   *      Cx >= d
   */
  Solution solveQP(const Matrix& Q,
                   const Vector& c,
                   const Matrix& A,
                   const Vector& b,
                   const Matrix& C,
                   const Vector& d);
}

#endif //QP_H
