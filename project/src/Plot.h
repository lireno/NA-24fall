#ifndef PLOT_H
#define PLOT_H

#include <fstream>
#include <iostream>
#include <string>

template <typename T>
void savePlottingData(const T& f, const std::string& filename, size_t numSample = 300) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    f.plot(outFile, numSample);
    outFile.close();
}

template <typename T>
void savePlottingData(const T& f, const std::string& filename, double a, double b, size_t numSample = 300) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    f.plot(outFile, a, b, numSample);
    outFile.close();
}

#endif // PLOT_H
