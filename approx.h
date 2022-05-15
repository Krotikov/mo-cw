#pragma once

#ifndef APPROX_H
#define APPROX_H

#include "func.h"

class Approx {
public:

  Approx(double eps) {
    this->f = Func();
    this->eps = eps;
    numIter = 0;
  }

  Approx(double eps, Func &f) {
    this->f = f;
    this->eps = eps;
    numIter = 0;
  }

  size_t numOfCalls() {
    return f.numOfCalls();
  }

  size_t getNumIter() {
    return numIter;
  }

  bool checkArgs(double x1, double x2, double x3, double f1, double f2, double f3) {
    bool t = false;

    if ((x1 < x2) && (x2 < x3))
      if ((f1 >= f2) && (f2 <= f3))
        t = true;

    return t;
  }

  double findSol() {
    double x1 = f.getL();
    double x3 = f.getR();
    double x2 = (x1 + x3) / 2;
    double res = -0.75;
    double f1 = f.calc(x1);
    double f2 = f.calc(x2);
    double f3 = f.calc(x3);

    do {
      //check
      bool tmp = checkArgs(x1, x2, x3, f1, f2, f3);
      double t = x2;
      double fc = f2;

      if (tmp == false) {
        //if (f3 < f1)
          //return x3;

        double xm1 = x1;
        double xm2 = x2;
        double xm3 = x3;
        double fm2 = f2;

        for (size_t p2 = 2; 2 * eps < (x3 - x1) / p2; p2 *= 2) {
          t = (xm1 + xm2) / 2;
          fc = f.calc(t);
          if (fc < fm2) {
            if (checkArgs(x1, t, x3, f1, fc, f3) == true)
              break;

            fm2 = fc;
            xm2 = t;
            xm3 = xm2;
            continue;
          }
          t = (xm2 + xm3) / 2;
          fc = f.calc(t);
          if (fc < fm2) {
            if (checkArgs(x1, t, x3, f1, fc, f3) == true) 
              break;
            
            fm2 = fc;
            xm2 = t;
            xm1 = xm2;
          }
        }
      }
      //assert(t);
      x2 = t;
      f2 = fc;

      double a1 = (f2 - f1) / (x2 - x1);
      double a2 = ((f3 - f1) / (x3 - x1) - a1) / (x3 - x2);

      double sol = (x1 + x2 - a1 / a2) / 2;

      double delta = abs(sol - res);
      res = sol;
      numIter++;

      if (delta < eps)
        break;

      if (sol > x2) {
        x1 = x2;
        x2 = sol;
        f1 = f2;
        f2 = f.calc(x2);
      }
      else {
        x3 = x2;
        x2 = sol;
        f3 = f2;
        f2 = f.calc(x2);
      }

    } while (1);

    return res;
  }


private:
  size_t numIter;
  double eps;
  Func f;
};

#endif //APPROX_H
