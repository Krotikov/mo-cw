#ifndef BISEC_H
#define BISEC_H

#include <limits.h>
#include <cfloat>
#include "func.h"

typedef struct {
  double a, b;
  double ans;
}ans_t;

double UniformSearchMethod(double a, double b, int n, double eps, Func& f) {
  double h = (b - a) / n;
  double xi = a;
  double min = DBL_MAX, minxa = a, minxb = b;

  double fa = f.calc(a), fb = f.calc(b), fprev = fa;
  ans_t ans;
  while (b - a > eps) {

    double y;

    for (int j = 0; j <= n; ++j, xi += h) {

      if (xi == a) 
        y = fa;
      else if (xi == b) 
        y = fb;
      else
        y = f.calc(xi);

      if (min >= y) {
        min = y;
        if (xi - h < a)
          minxa = a;
        else
          minxa = xi - h;

        if (xi + h > b)
          minxb = b;
        else
          minxb = xi + h;
      }
      //min < y->we can stop
      else {
        fb = y;
        break;
      }

      fa = fprev;
      fprev = y;
    }

    a = minxa;
    b = minxb;
    h = (b - a) / n;
    min = DBL_MAX;
    xi = a;

  }
  ans.a = a;
  ans.b = b;
  ans.ans = (b + a) / 2.0;

  f.numOfCalls();
  return ans.ans;
  //return (b + a) / 2.0;
}


#endif //BISEC_H