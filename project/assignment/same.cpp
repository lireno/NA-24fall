// This file is to test that PPSpline and BSpline has the same result

#include "BSpline.h"
#include "CurveFitting.h"
#include "PPSpline.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {

    std::vector<double> nodes = {1, 2, 4, 5, 7};
    std::vector<double> values = {6, 7, 6, 5, 6};

    NaturalCubicBSpline NaturalCubicBSpline(nodes, values);
    NaturalCubicPPSpline NaturalCubicPPSpline(nodes, values);

    PeriodicCubicBSpline PeriodicCubicBSpline(nodes, values);
    PeriodicCubicPPSpline PeriodicCubicPPSpline(nodes, values);

    CompleteCubicBSpline CompleteCubicBSpline(nodes, values, 2, 1);
    CompleteCubicPPSpline CompleteCubicPPSpline(nodes, values, 2, 1);

    std::cout << "x\tN_B\tN_PP\t"
              << "P_B\tP_PP\t"
              << "C_B\tC_PP" << std::endl;
    for (double x = 1; x < 6; x += 0.6) {
        std::cout << x << "\t"
                  << NaturalCubicBSpline(x) << "\t"
                  << NaturalCubicPPSpline(x) << "\t"
                  << PeriodicCubicBSpline(x) << "\t"
                  << PeriodicCubicPPSpline(x) << "\t"
                  << CompleteCubicBSpline(x) << "\t"
                  << CompleteCubicPPSpline(x) << std::endl;
    }

    return 0;
}
