#include "BSpline.h"
#include "CurveFitting.h"
#include "Function.h"
#include "Plot.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class c1 : public Curve {
  public:
    c1() : Curve(0, 2 * Pi) {}
    Point operator()(double t) const override {
        double x = sqrt(3) * cos(t);
        double y = 2.0 / 3.0 * (sqrt(std::abs(x)) + sqrt(3) * sin(t));
        return Point(x, y);
    }
};

class c2 : public Curve {
  public:
    c2() : Curve(0, 6 * Pi) {}
    Point operator()(double t) const override {
        return Point(sin(t) + t * cos(t), cos(t) - t * sin(t));
    }
};

int main() {
    c1 c1_curve;
    c2 c2_curve;

    // CurveFitting cf1_c1(c1_curve, 10);
    // savePlottingData(cf1_c1, "curve_fitting_c1_10.txt");

    // CurveFitting cf2_c1(c1_curve, 40);
    // savePlottingData(cf2_c1, "curve_fitting_c1_40.txt");

    // CurveFitting cf3_c1(c1_curve, 160);
    // savePlottingData(cf3_c1, "curve_fitting_c1_160.txt");

    CurveFitting eqs_cf1_c1(c2_curve, 1.0);
    savePlottingData(eqs_cf1_c1, "eqs_curve_fitting_c1_0.1.txt");
    std::cout << "Control points size: " << eqs_cf1_c1.control_points_size() << std::endl;

    CurveFitting eqs_cf2_c1(c2_curve, 3.0);
    savePlottingData(eqs_cf2_c1, "eqs_curve_fitting_c1_0.01.txt");
    std::cout << "Control points size: " << eqs_cf2_c1.control_points_size() << std::endl;

    // CurveFitting cf1_c2(c2_curve, 10);
    // savePlottingData(cf1_c2, "curve_fitting_c2_10.txt");

    // CurveFitting cf2_c2(c2_curve, 40);
    // savePlottingData(cf2_c2, "curve_fitting_c2_40.txt");

    // CurveFitting cf3_c2(c2_curve, 160);
    // savePlottingData(cf3_c2, "curve_fitting_c2_160.txt");

    return 0;
}