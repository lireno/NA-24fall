#include "BSpline.h"
#include "PPSpline.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void saveSplineData(const BSpline& spline, const std::string& filename, size_t numSamplesPerInterval = 100) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    spline.plotSpline(outFile, numSamplesPerInterval);
    outFile.close();
    std::cout << "Saved spline points to " << filename << std::endl;
}

int main() {
    // Define nodes and values for multiple splines
    std::vector<std::vector<double>> nodes = {
        {0, 1, 2, 3, 4, 5},
        {0, 1, 2, 3, 4, 5},
        {0, 1, 2, 3, 4, 5}};
    std::vector<std::vector<double>> values = {
        {0, 1, 4, 9, 16, 25},
        {0, -1, -4, -9, -16, -25},
        {0, 2, 8, 18, 32, 50}};

    // Output directory for plot data
    std::string outputDir = "./plotting/data/";

    // Create splines and save data
    for (size_t i = 0; i < nodes.size(); ++i) {
        // Create a natural cubic spline
        NaturalCubicBSpline spline(nodes[i], values[i]);

        // Save spline data to a file
        std::string filename = outputDir + "spline_curve_" + std::to_string(i) + ".txt";
        saveSplineData(spline, filename);
    }

    return 0;
}
