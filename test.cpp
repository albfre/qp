#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "qp.h"
#include "MathUtil.h"

void print(const MathUtil::Vector& x) {
  std::stringstream ss;
  for (const auto& xi : x) {
    ss << xi << ", ";
  }
  ss << std::endl;
  std::cout << ss.str();
}

int main( int argc, const char* argv[] )
{
  const size_t xSize = 100;
  const size_t ySize = 20;
  MathUtil::Matrix H(xSize, MathUtil::Vector(xSize));
  MathUtil::Vector g(xSize, 1);
  for (size_t i = 0; i < xSize; ++i) {
    H[i][i] = 1.0;
  }

  MathUtil::Matrix A(ySize, QP::Vector(xSize));
  MathUtil::Vector b(ySize);
  for (size_t i = 0; i < ySize; ++i) {
    A[i][i] = -2.0;
    A[i][i + 1] = -1.0;
    b[i] = 0.2;
  }
  const auto s = QP::solveQP(H, g, A, b);
  print(s.x);
  std::cout << "Objective value: " << s.objectiveValue << std::endl;
  return 0;
}
