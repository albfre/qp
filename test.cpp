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


int main(int argc, const char* argv[])
{
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
  return 0;
}
