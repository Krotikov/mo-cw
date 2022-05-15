#ifndef GRAD_H
#define GRAD_H

#include "func2.h"
#include "func.h"
#include "approx.h"
#include "bisection.h"
#include <iostream>


class GradMet {
public:

  GradMet(double eps, double r) {
    this->eps = eps;
    this->f = Func2(r);
    res = new double[2];
    res[0] = 0;
    res[1] = 0;
    n = 0;
    g = 0;
  }

  GradMet(double eps, Func2 f) {
    this->eps = eps;
    this->f = f;
    res = new double[2];
    res[0] = res[1] = 0;
    n = 0;
    g = 0;
  }

  size_t numOfCalls() {
    return n;
  }

  size_t numOfGrad() {
    return g;
  }

  double maxAbs(double x1, double x2) {
    if (abs(x1) > abs(x2))
      return abs(x1);

    return abs(x2);
  }

  double* calc() {
    double* xk = res;
    double* pr = new double[2];
    size_t k = 1;
    double gradNorm = 1;
    double gPrev = 1;
    double dgrad = 1;
    double fine = 1;
    double p = 1;
    double r = f.getR();

    double x1 = 1 / sqrt(2);
    double x2 = 1 / sqrt(2);

    do {
      do {
        double* d = f.grad(xk[0], xk[1]);
        r = f.getR();
        Func func(1e-12, 1, xk[0], xk[1], d[0], d[1], r);

        double a = UniformSearchMethod(1e-12, 1, 8, eps*0.1, func);
        n += func.numOfCalls();
        if (a == -1)
          break;

        pr[0] = xk[0];
        pr[1] = xk[1];
        xk[0] = xk[0] - a * d[0];
        xk[1] = xk[1] - a * d[1];

        gPrev = gradNorm;
        gradNorm = f.gradNorm(xk[0], xk[1]);
        fine = f.fine(xk[0], xk[1], eps);
        dgrad = (gradNorm - gPrev) * (gradNorm - gPrev);

        std::cout << k << ") " << "Ans x1 = " << xk[0] << ", x2 = " << xk[1] << ", grad = " << gradNorm
          << ", qk = " << maxAbs(xk[0] - x1, xk[1] - x2) / pow(maxAbs(pr[0] - x1, pr[1] - x2), p)
          << ", ak = " << a << ", fine = " << fine << ", dgrad = " << dgrad << std::endl;
        ++k;

        if (k > 110) {
          std::cout << "CAPUT, GG WELL PLAYED" << std::endl;
          break;
        }

      } while (dgrad > eps);
    } while (fine > eps);

    n += f.numOfCalls();
    g += f.numOfGrad();
    res = xk;

    return res;
  }

  ~GradMet() {
    delete[] res;
  }


private:
  double* res;
  double eps;
  Func2 f;
  size_t n;
  size_t g;
};


#endif //GRAD_H
