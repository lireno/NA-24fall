#include "BSpline.h"
#include "Function.h"
#include "PPSpline.h"
#include "Plot.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Define the function f(x) = 1 / (1 + 25 * x^2)
class F1 : public Function {
  public:
    double operator()(double x) const override {
        return 1.0 / (1.0 + x * x);
    }
};

// Main function
int main() {
    std::vector<double> nodes = {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    F1 f;
    CompleteCubicBSpline s1(f, nodes);
    quarticBSpline s2(f, -5, 5);
    std::string outputDir = "./plotting/data/";
    savePlottingData(s1, "completeCubicBSpline.txt", 10);
    savePlottingData(s2, "quarticBSpline.txt");
    savePlottingData(f, "f1.txt", -5.0, 5.0);
    return 0;
}
