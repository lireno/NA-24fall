#include "BSpline.h"
#include "CurveFitting.h"
#include "PPSpline.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {

    std::vector<double> nodes = {1, 2, 4, 5, 7};
    std::vector<double> a = {2, 4, 6, 7, 6, 5, 6, 2};

    BSpline bspline(nodes, a, 4);

    for (double x = 1; x < 6; x += 0.6) {
        std::cout << x << "\t" << bspline(x) << std::endl;
    }

    return 0;
}