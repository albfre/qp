#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "qp.h"

int main( int argc, const char* argv[] )
{
  const size_t xSize = 10;
  const size_t ySize = 1;
  QP::Matrix H(xSize, QP::Vector(xSize));
  QP::Vector g(xSize);
  for (size_t i = 0; i < xSize; ++i) {
    H[i][i] = 1.0;
  }

  QP::Matrix A(ySize, QP::Vector(xSize));
  QP::Vector b(ySize);
  for (size_t i = 0; i < ySize; ++i) {
    A[i][i] = -2.0;
    A[i][i + 1] = -1.0;
    b[i] = 0.2;
  }
  const auto s = QP::solveQP(H, g, A, b);
  std::cout << "Objective value: " << s.objectiveValue << std::endl;
  for (const auto& xi : s.x) {
    std::cout << xi << ", ";
  }
  std::cout << std::endl;
  return 0;
}
