#include "Function.hpp"
#include "Interpolator.hpp"
#include <cmath>
#include <iostream>
#include <vector>

class RungeFunction : public Function {
  public:
    double operator()(double x) const override {
        return 1.0 / (1.0 + x * x);
    }
};

int main() {
    double a = -5.0;
    double b = 5.0;
    std::vector<int> n_values = {2, 4, 6, 8};

    RungeFunction runge;

    for (int n : n_values) {
        std::vector<double> x_values;
        for (int i = 0; i <= n; ++i) {
            double x_i = a + (b - a) * i / n;
            x_values.push_back(x_i);
        }

        NewtonInterpolator interpolator(runge, x_values);

        std::cout << "The main dialog of divided diffenrences table and x_i for n =" << n << ":\n";
        interpolator.printInternalData();
        std::cout << "\n";
    }
    return 0;
}