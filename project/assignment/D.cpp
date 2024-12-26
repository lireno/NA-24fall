#include "BSpline.h"
#include "Function.h"
#include "PPSpline.h"
#include "Plot.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Define the function f(x) = 1 / (1 + x^2)
class F1 : public Function {
  public:
    double operator()(double x) const override {
        return 1.0 / (1.0 + x * x);
    }
};

void reportInterpolationErrors(const Function& f, const Function& spline, const std::vector<double>& testSites) {
    std::cout << std::setw(10) << "x"
              << std::setw(20) << "f(x)"
              << std::setw(20) << "S(x)"
              << std::setw(20) << "|S(x) - f(x)|" << std::endl;

    for (double x : testSites) {
        double fx = f(x);
        double sx = spline(x);
        double error = std::abs(sx - fx);

        std::cout << std::setw(10) << x
                  << std::setw(20) << fx
                  << std::setw(20) << sx
                  << std::setw(20) << error << std::endl;
    }
}

int main() {
    // Define nodes and function
    std::vector<double> nodes = {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    F1 f;

    // Create spline interpolants
    CompleteCubicBSpline s1(f, nodes); // Complete cubic B-spline
    quarticBSpline s2(f, -5, 5);       // Quartic B-spline

    // Define test sites
    std::vector<double> testSites = {-3.5, -3.0, -0.5, 0.0, 0.5, 3.0, 3.5};

    // Report interpolation errors
    std::cout << "Errors for Complete Cubic B-Spline:\n";
    reportInterpolationErrors(f, s1, testSites);

    std::cout << "\nErrors for Quartic B-Spline:\n";
    reportInterpolationErrors(f, s2, testSites);

    return 0;
}
