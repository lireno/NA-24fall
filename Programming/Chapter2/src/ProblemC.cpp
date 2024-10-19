#include "Function.hpp"
#include "Interpolator.hpp"
#include <cmath>
#include <iostream>
#include <vector>

std::vector<double> generateChebyshevNodes(int n) {
    std::vector<double> nodes(n);
    for (int i = 0; i < n; ++i) {
        nodes[i] = cos((2.0 * i + 1.0) / (2.0 * n) * M_PI);
    }
    return nodes;
}

class RungeFunction : public Function {
  public:
    double operator()(double x) const override {
        return 1.0 / (1.0 + 25 * x * x);
    }
};

int main() {
    std::vector<int> n_values = {5, 10, 15, 20};

    RungeFunction runge;

    for (int n : n_values) {
        std::vector<double> chebyshevNodes = generateChebyshevNodes(n);
        NewtonInterpolator interpolator(runge, chebyshevNodes);

        std::cout << "The main dialog of divided diffenrences table and x_i for n =" << n << ":\n";
        interpolator.printInternalData();
        std::cout << "\n";
    }
    return 0;
}