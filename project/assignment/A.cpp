#include "Function.h"
#include "PPSpline.h"
#include "Plot.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Define the function f(x) = 1 / (1 + 25 * x^2)
class RungeFunction : public Function {
  public:
    double operator()(double x) const override {
        return 1.0 / (1.0 + 25.0 * x * x);
    }
};

// Compute the maximum error at midpoints
double computeMaxError(const NaturalCubicPPSpline& spline, const std::vector<double>& nodes, const Function& f) {
    double maxError = 0.0;

    for (size_t i = 0; i < nodes.size() - 1; ++i) {
        double midpoint = (nodes[i] + nodes[i + 1]) / 2.0;
        double error = std::abs(spline(midpoint) - f(midpoint));
        maxError = std::max(maxError, error);
    }

    return maxError;
}

// Generate evenly spaced nodes in the interval [a, b]
std::vector<double> generateNodes(double a, double b, int N) {
    std::vector<double> nodes(N);
    double step = (b - a) / (N - 1);
    for (int i = 0; i < N; ++i) {
        nodes[i] = a + i * step;
    }
    return nodes;
}

int main() {
    RungeFunction f;
    double a = -1.0, b = 1.0;
    std::string outputDir = "./plotting/data/";
    std::vector<int> N_values = {6, 11, 21, 41, 81};

    std::vector<double> maxErrors;
    for (int N : N_values) {
        std::vector<double> nodes = generateNodes(a, b, N);
        std::vector<double> values;
        for (double x : nodes) {
            values.push_back(f(x));
        }

        NaturalCubicPPSpline spline(nodes, values);

        double maxError = computeMaxError(spline, nodes, f);
        maxErrors.push_back(maxError);

        std::string filename = "splines with N = " + std::to_string(N) + ".txt";
        savePlottingData(spline, filename);
    }

    std::cout << std::setw(10) << "N"
              << std::setw(20) << "Max Error" << std::endl;

    for (int i = 0; i <= 4; ++i) {

        std::cout << std::setw(10) << N_values[i]
                  << std::setw(20) << maxErrors[i] << std::endl;
    }

    savePlottingData(f, "runge_function.txt", 1.0, -1.0);
    return 0;
}
