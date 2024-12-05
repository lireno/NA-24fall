#include "BSpline.h"
#include "CurveFitting.h"
#include "Function.h"
#include "Plot.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class spcurve : public SphereCurve {
  public:
    spcurve() : SphereCurve(0, 2 * Pi) {}
    Point3D operator()(double t) const override {
        double u = std::cos(t);
        double v = std::sin(t);
        return Point3D(std::sin(u) * std::cos(v), std::sin(u) * std::sin(v), std::cos(u));
    }
};

int main() {
    spcurve spcurve;
    SphereCurvefit cf1(spcurve, {0, 0, -1}, 10);
    SphereCurvefit cf2(spcurve, {0, 0, -1}, 50);
    SphereCurvefit cf3(spcurve, {0, 0, -1}, 100);

    savePlottingData(cf1, "sphere_curve_fitting_10.txt");
    savePlottingData(cf2, "sphere_curve_fitting_50.txt");
    savePlottingData(cf3, "sphere_curve_fitting_100.txt");
    savePlottingData(spcurve, "sphere_curve.txt");
    return 0;
}