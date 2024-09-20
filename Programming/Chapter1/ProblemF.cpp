#include "EquationSolver.hpp"
#include "Function.hpp"
#include <cmath>
#include <iostream>

const double Pi = acos(-1.);

class NoseInFailure : public Function {
  private:
    double A, B, C, E;

  public:
    NoseInFailure(double l, double h, double D, double beta1) {
        A = l * sin(beta1);
        B = l * cos(beta1);
        C = (h + 0.5 * D) * sin(beta1) - 0.5 * D * tan(beta1);
        E = (h + 0.5 * D) * cos(beta1) - 0.5 * D;
    }

    double operator()(double alpha) const {
        return A * sin(alpha) * cos(alpha) + B * sin(alpha) * sin(alpha) - C * cos(alpha) - E * sin(alpha);
    }
};

void solve_a() {
    double l = 89.0, h = 49.0, D = 55.0, beta1 = 11.5 * (M_PI / 180.0);
    std::cout << "Finding the Maximum with D = 55 in" << std::endl;
    NoseInFailure failure_eq(l, h, D, beta1);
    Newton_Method solver(failure_eq, 33.0 * (M_PI / 180.0));
    double alpha = solver.solve();
    std::cout << "Maximum angle alpha: " << alpha * (180.0 / M_PI) << " degrees." << std::endl;
}

void solve_bc() {
    double l = 89.0, h = 49.0, D = 30.0, beta1 = 11.5 * (M_PI / 180.0);
    std::cout << "Finding the Maximum with D = 30 in" << std::endl;
    NoseInFailure failure_eq(l, h, D, beta1);
    Newton_Method solver(failure_eq, 33.0 * (M_PI / 180.0));
    double alpha = solver.solve();
    std::cout << "Maximum angle alpha: " << alpha * (180.0 / M_PI) << " degrees." << std::endl;

    std::cout << "Finding the Maximum with D = 30 and initial values of 75 and 100 degree" << std::endl;
    Secant_Method solver2(failure_eq, 75.0 * (M_PI / 180.0), 100.0 * (M_PI / 180.0));
    double alpha2 = solver2.solve();
    std::cout << "Maximum angle alpha: " << alpha2 * (180.0 / M_PI) << " degrees." << std::endl;

    std::cout << "Finding the Maximum with D = 30 and initial values of 130 and 140 degree" << std::endl;
    Secant_Method solver3(failure_eq, 130.0 * (M_PI / 180.0), 140.0 * (M_PI / 180.0));
    double alpha3 = solver3.solve();
    std::cout << "Maximum angle alpha: " << alpha3 * (180.0 / M_PI) << " degrees." << std::endl;

    std::cout << "Finding the Maximum with D = 30 and initial values of 160 and 175 degree" << std::endl;
    Secant_Method solver4(failure_eq, 160.0 * (M_PI / 180.0), 175.0 * (M_PI / 180.0));
    double alpha4 = solver4.solve();
    std::cout << "Maximum angle alpha: " << alpha4 * (180.0 / M_PI) << " degrees." << std::endl;
}

int main() {
    solve_a();
    solve_bc();
    return 0;
}