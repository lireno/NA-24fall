#include "Interpolator.hpp"
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    std::vector<double> time = {0, 3, 5, 8, 13};
    std::vector<double> displacement = {0, 225, 383, 623, 993};
    std::vector<double> velocity = {75, 77, 80, 74, 72};

    HermiteInterpolator interpolator(time, displacement, velocity);

    // std::cout << "the divided differences are following:" << std::endl;
    // interpolator.printInternalData();

    double t_target = 10;
    double displacement_10 = interpolator.interpolate(t_target);
    std::cout << "(a) The position at t = 10 seconds: " << displacement_10 << " feet" << std::endl;

    double step = 0.001;
    double old_displacement = 225.0;
    double max_velocity = 80.0;
    double t_max_velocity = 5.0;
    for (double t = 3.001; t <= 8; t += step) { // clearly that the max velocity is between 3 and 8
        double new_displacement = interpolator.interpolate(t);
        double velocity = (new_displacement - old_displacement) / step;
        if (velocity > max_velocity) {
            max_velocity = velocity;
            t_max_velocity = t;
        }
        old_displacement = new_displacement;
    }

    std::cout << "(b) The maximum velocity is " << max_velocity << " feet/second, reached at t = " << t_max_velocity << " seconds" << std::endl;

    if (max_velocity > 81.0) {
        std::cout << " The maximum velocity exceeds 81 feet/second" << std::endl;
    } else {
        std::cout << " The maximum velocity does not exceed 81 feet/second" << std::endl;
    }

    return 0;
}