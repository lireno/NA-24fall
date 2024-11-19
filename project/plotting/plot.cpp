#include "BSpline.h"
#include "CurveFitting.h"
#include "PPSpline.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <typename T>
void savePlottingData(const T& f, const std::string& filename, size_t numSample = 300) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    f.plot(outFile, numSample);
    outFile.close();
    std::cout << "Saved spline points to " << filename << std::endl;
}

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
