#ifndef FUNCTION
#define FUNCTION
#include <cmath>

class Function {
  public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const {
        double h = 1e-6;
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h);
    }

    bool is_Continuous(double x_0, double eps = 1e-2) const {
        double h = 1e-6;
        double leftLimit = (*this)(x_0 - h);
        double rightLimit = (*this)(x_0 + h);
        double f_x = (*this)(x_0);
        return std::abs(leftLimit - f_x) < eps && std::abs(rightLimit - f_x) < eps;
    }
};

#endif