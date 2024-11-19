#ifndef PLOT_H
#define PLOT_H

#include <fstream>
#include <iostream>
#include <string>

// Template function to save plotting data for a function or spline
template <typename T>
void savePlottingData(const T& f, const std::string& filename, size_t numSample = 300) {
    std::string fullPath = "plotting/data/" + filename;
    std::ofstream outFile(fullPath);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << fullPath << " for writing!" << std::endl;
        return;
    }

    f.plot(outFile, numSample);
    outFile.close();

    std::cout << "Saved spline points to " << fullPath << std::endl;
}

template <typename T>
void savePlottingData(const T& f, const std::string& filename, double a, double b, size_t numSample = 300) {
    std::string fullPath = "plotting/data/" + filename;
    std::ofstream outFile(fullPath);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << fullPath << " for writing!" << std::endl;
        return;
    }

    f.plot(outFile, a, b, numSample);
    outFile.close();

    std::cout << "Saved spline points to " << fullPath << std::endl;
}

#endif // PLOT_H
