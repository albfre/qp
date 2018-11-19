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

int main(int argc, const char* argv[])
{
  const size_t xSize = 1000;
  const size_t ySize = 200;
  MathUtil::Matrix Q(xSize, MathUtil::Vector(xSize));
  MathUtil::Vector c(xSize, 1);
  for (size_t i = 0; i < xSize; ++i) {
    Q[i][i] = 3 * xSize;
    c[i] = i;
  }
  for (size_t i = 0; i < xSize; ++i) {
    for (size_t j = i + 1 ; j < xSize; ++j) {
      Q[i][j] = Q[j][i] = i + j;
    }
  }
  MathUtil::Matrix A(0, QP::Vector(0));
  MathUtil::Vector b(0);

  MathUtil::Matrix C(ySize, QP::Vector(xSize));
  MathUtil::Vector d(ySize);
  for (size_t i = 0; i < ySize; ++i) {
    for (size_t j = 0; j < xSize; ++j) {
      C[i][j] = i - j;
    }
    d[i] = i;
  }
  const auto s = QP::solveQP(Q, c, A, b, C, d);
//  print(s.x);
  std::cout << "Objective value: " << s.objectiveValue << std::endl;
  std::cout << "Infeasibility: " << s.infeasibility << std::endl;
  return 0;
}
