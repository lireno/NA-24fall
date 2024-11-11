#ifndef FUNCTION_H
#define FUNCTION_H
#include <cmath>

class Function {
  public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const {
        double h = 1e-6;
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h);
    }
};
#endif // FUNCTION_H
