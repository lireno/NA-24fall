#include "BSpline.h"
#include "Function.h"
#include "PPSpline.h"
#include "Plot.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

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

    savePlottingData(s1, "completeCubicBSpline.txt");
    savePlottingData(s2, "quarticBSpline.txt");
    savePlottingData(f, "RungeFunction.txt", -5.0, 5.0);
    return 0;
}
