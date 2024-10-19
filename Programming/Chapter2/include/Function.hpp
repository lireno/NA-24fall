#ifndef FUNCTION
#define FUNCTION
#include <cmath>

class Function {
  public:
    virtual double operator()(double x) const = 0;
};

#endif