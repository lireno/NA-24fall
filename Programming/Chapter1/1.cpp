#include "EquationSolver.hpp"
#include "Function.hpp"
#include <cmath>
#include <iostream>

const double Pi = acos(-1.);

class F2 : public Function {
  public:
    double operator()(double x) const {
        return sin(x) + pow(cos(x), 2) * log((1 + sin(x)) / cos(x));
    }
};

class F1 : public Function {
  public:
    double operator()(double x) const {
        return pow(cos(x), 3) - (1 + sin(x)) * cos(x) * sin(x) * (2 * log(tan(x) + 1 / cos(x)) - 1) - 1;
    }
};

void solve_f1() {
    std::cout << "Solving x^{-1} - \\tan x on [0, \\pi/2]" << std::endl;
    Bisection_Method solver_f1(F1(), 0.01, 1.5);
    double x = solver_f1.solve();
    std::cout << "A root is: " << x << std::endl;
}

int main() {
    // Output values of F2 from 0 to Pi/2, taking 30 samples
    std::cout << "Values of F2 from 0 to Pi/2:" << std::endl;
    for (double i = 0; i < 180; i += 1) {
        std::cout << "f(" << (56 + i / 60) * Pi / 180 << ") = " << F2()((56 + i / 60) * Pi / 180) << "\t";
    }
    std::cout << std::endl;

    solve_f1();
    return 0;
}