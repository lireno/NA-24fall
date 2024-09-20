#include "EquationSolver.hpp"
#include "Function.hpp"
#include <cmath>
#include <iostream>

const double Pi = acos(-1.);

class F1 : public Function {
  public:
    double operator()(double x) const {
        return 1.0 / x - tan(x);
    }
};

class F2 : public Function {
  public:
    double operator()(double x) const {
        return 1.0 / x - pow(2.0, x);
    }
};

class F3 : public Function {
  public:
    double operator()(double x) const {
        return pow(2.0, -x) + exp(x) + 2 * cos(x) - 6;
    }
};

class F4 : public Function {
  public:
    double operator()(double x) const {
        return (((x + 4.0) * x + 3.0) * x + 5.0) / (((2 * x - 9.0) * x + 18.0) * x - 2.0);
    }
};

void solve_f1() {
    std::cout << "Solving x^{-1} - \\tan x on [0, \\pi/2]" << std::endl;
    Bisection_Method solver_f1(F1(), 0, Pi / 2);
    double x = solver_f1.solve();
    std::cout << "A root is: " << x << std::endl;
}

void solve_f2() {
    std::cout << "Solving x^2 - 2 on [0, 1]" << std::endl;
    Bisection_Method solver_f2(F2(), 0, 1);
    double x = solver_f2.solve();
    std::cout << "A root is: " << x << std::endl;
}

void solve_f3() {
    std::cout << "Solving 2^{-x} - e^x + 2 \\cos x - 6 on [1, 3]" << std::endl;
    Bisection_Method solver_f3(F3(), 1, 2);
    double x = solver_f3.solve();
    std::cout << "A root is: " << x << std::endl;
}

void solve_f4() {
    std::cout << "Solving (x^3 + 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2) on [0, 4]" << std::endl;
    Bisection_Method solver_f4(F4(), 0, 4);
    double x = solver_f4.solve();
    std::cout << "A root is: " << x << std::endl;
}

int main() {
    solve_f1();
    solve_f2();
    solve_f3();
    // solve_f4();
    return 0;
}