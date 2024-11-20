// This file is to test that PPSpline and BSpline has the same result

#include "BSpline.h"
#include "CurveFitting.h"
#include "PPSpline.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class F1 : public Function {
  public:
    virtual double operator()(double x) const override {
        return 1 / (1 + x * x);
    };
};

int main() {

    std::vector<double> nodes = {1, 2, 4, 5, 7};
    std::vector<double> values = {6, 7, 6, 5, 6};

    NaturalCubicBSpline naturalCubicBSpline(nodes, values);
    NaturalCubicPPSpline naturalCubicPPSpline(nodes, values);

    PeriodicCubicBSpline periodicCubicBSpline(nodes, values);
    PeriodicCubicPPSpline periodicCubicPPSpline(nodes, values);

    CompleteCubicBSpline completeCubicBSpline(nodes, values, 2, 1);
    CompleteCubicPPSpline completeCubicPPSpline(nodes, values, 2, 1);

    std::cout << std::setw(5) << "x" << "\t"
              << std::setw(10) << "N_B" << "\t"
              << std::setw(10) << "N_PP" << "\t"
              << std::setw(10) << "P_B" << "\t"
              << std::setw(10) << "P_PP" << "\t"
              << std::setw(10) << "C_B" << "\t"
              << std::setw(10) << "C_PP" << std::endl;

    for (double x = 1; x < 6; x += 0.6) {
        std::cout << std::setw(5) << x << "\t"
                  << std::setw(10) << naturalCubicBSpline(x) << "\t"
                  << std::setw(10) << naturalCubicPPSpline(x) << "\t"
                  << std::setw(10) << periodicCubicBSpline(x) << "\t"
                  << std::setw(10) << periodicCubicPPSpline(x) << "\t"
                  << std::setw(10) << completeCubicBSpline(x) << "\t"
                  << std::setw(10) << completeCubicPPSpline(x) << std::endl;
    }

    std::cout << " ------------------------------------------------------------------------------------------" << std::endl;

    F1 f;
    std::vector<double> nodes2 = {-5, -4, -3, -1, 0, 1, 3, 4, 5};

    NaturalCubicBSpline NaturalCubicBSpline2(f, nodes2);
    NaturalCubicPPSpline NaturalCubicPPSpline2(f, nodes2);

    PeriodicCubicBSpline PeriodicCubicBSpline2(f, nodes2);
    PeriodicCubicPPSpline PeriodicCubicPPSpline2(f, nodes2);

    CompleteCubicBSpline CompleteCubicBSpline2(f, nodes2);
    CompleteCubicPPSpline CompleteCubicPPSpline2(f, nodes2);

    std::cout << std::setw(5) << "x" << "\t"
              << std::setw(10) << "N_B" << "\t"
              << std::setw(10) << "N_PP" << "\t"
              << std::setw(10) << "P_B" << "\t"
              << std::setw(10) << "P_PP" << "\t"
              << std::setw(10) << "C_B" << "\t"
              << std::setw(10) << "C_PP" << std::endl;

    for (double x = -5; x < 5; x += 1.2) {
        std::cout << std::setw(5) << x << "\t"
                  << std::setw(10) << NaturalCubicBSpline2(x) << "\t"
                  << std::setw(10) << NaturalCubicPPSpline2(x) << "\t"
                  << std::setw(10) << PeriodicCubicBSpline2(x) << "\t"
                  << std::setw(10) << PeriodicCubicPPSpline2(x) << "\t"
                  << std::setw(10) << CompleteCubicBSpline2(x) << "\t"
                  << std::setw(10) << CompleteCubicPPSpline2(x) << std::endl;
    }
    return 0;
}
