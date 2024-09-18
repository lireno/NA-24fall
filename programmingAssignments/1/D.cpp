#include "EquationSolver.hpp"
#include "Function.hpp"
#include <cmath>
#include <iostream>

const double Pi = acos(-1.);

class F1 : public Function {
  public:
    double operator()(double x) const {
        return sin(x / 2) - 1;
    }
};

class F2 : public Function {
  public:
    double operator()(double x) const {
        return exp(x) - tan(x);
    }
};

class F3 : public Function {
  public:
    double operator()(double x) const {
        return ((x - 12) * x + 3) * x + 1;
    }
};

void solve1() {
    std::cout << "Finding the root of \\sin(x/2) - 1" << std::endl;
    Secant_Method solver(F1(), 0, Pi / 2);
    double x = solver.solve();
    std::cout << "A root is: " << x << std::endl;
}

void solve2() {
    std::cout << "Finding the root of e^x - \\tanx" << std::endl;
    Secant_Method solver(F2(), 1, 1.4);
    double x = solver.solve();
    std::cout << "A root is: " << x << std::endl;
}

void solve3() {
    std::cout << "Finding the root of x^3 - 12x^2 + 3x + 1 near 0" << std::endl;
    Secant_Method solver(F3(), 0, -0.5);
    double x = solver.solve();
    std::cout << "A root is: " << x << std::endl;
}

int main() {
    solve1();
    solve2();
    solve3();
    return 0;
}