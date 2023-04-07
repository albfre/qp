#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "QP.h"
#include "MathUtil.h"

void print(const MathUtil::Vector& x) {
  std::stringstream ss;
  for (const auto& xi : x) {
    ss << xi << ", ";
  }
  ss << std::endl;
  std::cout << ss.str();
}

void solveLP() {
  const size_t nx = 10;
  MathUtil::Matrix Q(nx, MathUtil::Vector(nx));
  MathUtil::Vector c(nx, 1);
  for (size_t i = 0; i < nx; ++i) {
    c[i] = i + 1;
  }
}

void testProblem() {
  const size_t xSize = 1000;
  const size_t ySize = 200;
  const auto print = false;
  MathUtil::Matrix Q(xSize, MathUtil::Vector(xSize));
  MathUtil::Vector c(xSize, 1);
  for (size_t i = 0; i < xSize; ++i) {
    for (size_t j = i + 1 ; j < xSize; ++j) {
      Q[i][j] = Q[j][i] = i + j;
    }
  }
  for (size_t i = 0; i < xSize; ++i) {
    Q[i][i] = 100 * xSize;
    c[i] = i;
  }

  if (print) {
    for (size_t i = 0; i < xSize; ++i) {
      for (size_t j = 0; j < xSize; ++j) {
        std::cout << Q[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }
  MathUtil::Matrix A(0, QP::Vector(0));
  MathUtil::Vector b(0);

  MathUtil::Matrix C(ySize, QP::Vector(xSize));
  MathUtil::Vector d(ySize);
  for (size_t i = 0; i < ySize; ++i) {
    for (size_t j = i; j < i + 4; ++j) {
      C[i][j] = static_cast<int>(i) - static_cast<int>(j);
    }
    d[i] = i;
  }
  if (print) {
    std::cout << "C:" << std::endl;
    for (size_t i = 0; i < ySize; ++i) {
      for (size_t j = 0; j < xSize; ++j) {
      std::cout << C[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }
  const auto printProgress = true;
  const auto s = QP::solveQP(Q, c, A, b, C, d, printProgress);
//  print(s.x);
  std::cout << "Objective value: " << std::setprecision(17) << s.objectiveValue << std::endl;
  std::cout << "Infeasibility: " << s.infeasibility << std::endl;
}

void testProblem2() {
  const size_t xSize = 2;
  const auto print = true;
  MathUtil::Matrix Q(xSize, MathUtil::Vector(xSize));
  MathUtil::Vector c(xSize, 1);
  const auto t0 = 1300.0;
  const auto a00 = 809.0;
  const auto a01 = 359.0;

  const auto t1 = 50.0;
  const auto a10 = 25.0;
  const auto a11 = 77.0;
  const auto denom = t0 * t0;

  Q[0][0] = (a00 * a00) / denom;
  Q[1][1] = (a01 * a01) / denom;
  Q[0][1] = (a00 * a01) / denom;
  Q[1][0] = (a00 * a01) / denom;
  c[0] = -t0 * a00 / denom;
  c[1] = -t0 * a01 / denom;


  if (print) {
    for (size_t i = 0; i < xSize; ++i) {
      for (size_t j = 0; j < xSize; ++j) {
        std::cout << Q[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }
  MathUtil::Matrix A(1, QP::Vector(xSize));
  MathUtil::Vector b(1);
  A[0][0] = 1;
  A[0][1] = 1;
  b[0] = 1;

  MathUtil::Matrix C(xSize, QP::Vector(xSize));
  MathUtil::Vector d(xSize);
  C[0][0] = 1;
  C[0][1] = 1;

  const auto printProgress = true;
  const auto s = QP::solveQP(Q, c, A, b, C, d, printProgress);
//  print(s.x);
  std::cout << "Objective value: " << std::setprecision(17) << s.objectiveValue << std::endl;
  std::cout << "Infeasibility: " << s.infeasibility << std::endl;
}

int main(int argc, const char* argv[])
{
  testProblem2();
  return 0;
}
