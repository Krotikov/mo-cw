#pragma once

#ifndef FUNC2_H
#define FUNC2_H

#include <cmath>

class Func2 {
public:
  Func2() {
    n = 0;
    g = 0;
    r = 1e+100;
    d = new double[2];
    d[0] = d[1] = 0;
  }

  Func2(double r) {
    this->r = r;
    n = 0;
    g = 0;
    d = new double[2];
    d[0] = d[1] = 0;
  }

  double g1(double x1, double x2) {
    double res = -x1;
    return (res < 0) ? 0 : res;
  }

  double derg1(double x1, double x2) {
    if (x1 < 0)
      return -1;

    return 0;
  }

  double g2(double x1, double x2) {
    double res = -x2;
    return (res < 0) ? 0 : res;
  }

  double derg2(double x1, double x2) {
    if (x2 < 0)
      return -1;

    return 0;
  }

  double g3(double x1, double x2) {
    double res = x1 * x1 + x2 * x2 - 3;
    return (res < 0) ? 0 : res;
  }

  double derg3x1(double x1, double x2) {
    if (x1 * x1 + x2 * x2 > 3)
      return 2 * x1;

    return 0;
  }

  double derg3x2(double x1, double x2) {
    if (x1 * x1 + x2 * x2 > 3)
      return 2 * x2;

    return 0;
  }

  double calc(double x1, double x2) {
    double res = (x1 - 1) * (x1 - 1) + (x2 - 1) * (x2 - 1) + 1/r*(g1(x1, x2) + g2(x1, x2) + g3(x1, x2));
    ++n;
    return res;
  }

  double fine(double x1, double x2, double eps) {
    double res = 1 / r * (g1(x1, x2) + g2(x1, x2) + g3(x1, x2));
    if(res > eps)
      r = r / sqrt(1.01);

    return res;
  }

  double* grad(double x1, double x2) {
    d[0] = 2 * (x1 - 1) + 1 / r * (derg1(x1, x2) + derg3x1(x1, x2));
    d[1] = 2 * (x2 - 1) + 1 / r * (derg2(x1, x2) + derg3x2(x1, x2));
    ++g;

    return d;
  }

  double gradNorm(double x1, double x2) {
    if (abs(d[0]) > abs(d[1]))
      return abs(d[0]);

    return abs(d[1]);
  }

  double* direction(double x1, double x2) {
    double h11 = 4 + 2 * exp(x1 * x1 + x2 * x2) + 4 * x1 * x1 * exp(x1 * x1 + x2 * x2);
    double h12 = 4 * x1 * x2 * exp(x1 * x1 + x2 * x2);
    double h21 = h12;
    double h22 = 2 + 2 * exp(x1 * x1 + x2 * x2) + 4 * x2 * x2 * exp(x1 * x1 + x2 * x2);
    double det = h11 * h22 - h12 * h21;

    d = grad(x1, x2);

    double* p = new double[2];
    p[0] = h22 / det * d[0] - h12 / det * d[1];
    p[1] = -h21 / det * d[0] + h11 / det * d[1];

    return p;
  }

  double getR() {
    return r;
  }

  size_t numOfCalls() {
    return n;
  }

  size_t numOfGrad() {
    return g;
  }

private:
  size_t n;
  size_t g;
  double r;
  double* d;
};

#endif //FUNC2_H