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
  MathUtil::Matrix H(xSize, MathUtil::Vector(xSize));
  MathUtil::Vector g(xSize, 1);
  for (size_t i = 0; i < xSize; ++i) {
    H[i][i] = 3 * xSize;
    g[i] = i;
  }
  for (size_t i = 0; i < xSize; ++i) {
    for (size_t j = i + 1 ; j < xSize; ++j) {
      H[i][j] = H[j][i] = i + j;
    }
  }

  MathUtil::Matrix A(ySize, QP::Vector(xSize));
  MathUtil::Vector b(ySize);
  for (size_t i = 0; i < ySize; ++i) {
    for (size_t j = 0; j < xSize; ++j) {
      A[i][j] = i - j;
    }
    b[i] = i;
  }
  const auto s = QP::solveQP(H, g, A, b);
  print(s.x);
  std::cout << "Objective value: " << s.objectiveValue << std::endl;
  std::cout << "Residual: " << s.residual << std::endl;
  return 0;
}
