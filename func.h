#pragma once

#ifndef FUNC_H
#define FUNC_H

#include <cmath>

class Func {
public:
  Func() {
    n = 0;
    r = 0;
    l = 0;
    c1 = c2 = d1 = d2 = 0;
    p = 10;
  }

  Func(double l, double r, double c1, double c2, double d1, double d2, double p) : 
    l(l), r(r), c1(c1), c2(c2), d1(d1), d2(d2), p(p) {
    n = 0;
  }

  double g1(double x1, double x2) {
    double res = -x1;
    return (res < 0) ? 0 : res;
  }

  double g2(double x1, double x2) {
    double res = -x2;
    return (res < 0) ? 0 : res;
  }

  double g3(double x1, double x2) {
    double res = x1 * x1 + x2 * x2 - 3;
    return (res < 0) ? 0 : res;
  }

  double calc(double x) {
    double x1 = c1 - x * d1;
    double x2 = c2 - x * d2;

    double res = (x1 - 1) * (x1 - 1) + (x2 - 1) * (x2 - 1) + 
                  1 / p * (g1(x1, x2) + g2(x1, x2) + g3(x1, x2)); //some func
    
    ++n;
    return res;
  }

  size_t numOfCalls() {
    return n;
  }

  double& getL() {
    return l;
  }

  double& getR() {
    return r;
  }

private:
  size_t n;
  double l;
  double r;
  double c1;
  double c2;
  double d1;
  double d2;
  double p;
};


#endif //FUNC_H