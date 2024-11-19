#ifndef FUNCTION_H
#define FUNCTION_H
#include <cmath>
#include <iostream>

class Function {
  public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const {
        double h = 1e-6;
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h);
    }

    void plot(std::ostream& os = std::cout, double xStart = 0, double xEnd = 1, size_t numSamples = 300) const {
        double step = (xEnd - xStart) / numSamples;
        for (size_t i = 0; i <= numSamples; ++i) {
            double x = xStart + i * step;
            os << x << " " << (*this)(x) << std::endl;
        }
    }
};

#endif // FUNCTION_H
