#include "Function.hpp"
#include "Interpolator.hpp"
#include <cmath>
#include <iostream>
#include <vector>

void predictSurvival(const std::vector<double>& days, const std::vector<double>& weights, const std::string& species, double day) {
    NewtonInterpolator interpolator(days, weights);
    // interpolator.printInternalData(); // this is for plot purposes
    double prediction = interpolator.interpolate(day);

    std::cout << "Prediction for " << species << " at day " << day << ": " << prediction << "\n";

    if (prediction <= 0) {
        std::cout << species << " larvae will likely die after " << day << " days." << std::endl;
    } else {
        std::cout << species << " larvae will likely survive after " << day << " days." << std::endl;
    }
}

int main() {
    std::vector<double> days = {0, 6, 10, 13, 17, 20, 28};
    std::vector<double> sp1_weights = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    std::vector<double> sp2_weights = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};

    double day_43 = 43;

    predictSurvival(days, sp1_weights, "Sp1", 43);
    predictSurvival(days, sp2_weights, "Sp2", 43);

    return 0;
}
