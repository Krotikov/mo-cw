#include "grad.h"
#include <iostream>


void Print(double* arr, size_t size) {
  for (size_t i = 0; i < size; ++i)
    std::cout << arr[i] << ", ";

  std::cout << std::endl;
}

void Print(size_t* arr, size_t size) {
  for (size_t i = 0; i < size; ++i)
    std::cout << arr[i] << ", ";

  std::cout << std::endl;
}


int main() {

  double eps[] = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
  double r[] = { 1e-2 };
  double x = 1;

  //Func2 f();
  std::cout << "x1 = x2 = " << x << std::endl;
  size_t size = sizeof(eps) / sizeof(double);
  size_t sizeR = sizeof(r) / sizeof(double);

  double* d = new double[size];
  size_t* nf = new size_t[size];
  size_t* ng = new size_t[size];
  size_t* sum = new size_t[size];

  for (size_t j = 0; j < sizeR; ++j) {
    for (size_t i = 0; i < size; ++i) {
      std::cout << "eps = " << eps[i] << std::endl;
      double* res;
      size_t n = 0;
      GradMet g(eps[i], r[j]);

      res = g.calc();

      d[i] = g.maxAbs(res[0] - x, res[1] - x);
      nf[i] = g.numOfCalls();
      ng[i] = g.numOfGrad();
      sum[i] = ng[i] + nf[i];

      //std::cout << "x1: " << res[0] << ", x2: " << res[1] << ", n: " << n << std::endl;
      std::cout << "================================" << std::endl;
    }
  }

  std::cout << "EPS: ";
  Print(eps, size);
  std::cout << "D: ";
  Print(d, size);
  std::cout << "NF: ";
  Print(nf, size);
  std::cout << "NG: ";
  Print(ng, size);
  std::cout << "SUM: ";
  Print(sum, size);

  std::cout << "================================" << std::endl;

  delete[] d;
  delete[] nf;
  delete[] ng;
  delete[] sum;

  return 0;
}