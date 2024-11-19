#include "Plot.h"
#include "BSpline.h"
#include "CurveFitting.h"
#include "PPSpline.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class C1 : public Curve {
  public:
    C1() : Curve(0, 2 * Pi) {}
    Point operator()(double t) const override {
        return {std::cos(t), std::sin(t)};
    }
};

int main() {

    std::string outputDir = "./plotting/data/";
    C1 c1;
    CurveFitting cf(c1, 0.1);
    savePlottingData(cf, outputDir + "curve_fitting.txt");

    return 0;
}
