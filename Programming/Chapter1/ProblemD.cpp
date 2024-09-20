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
    std::cout << "Finding the root of \\sin(x/2) - 1 with initial values of 0 and Pi/2" << std::endl;
    Secant_Method solver(F1(), 0, Pi / 2);
    double x = solver.solve();
    std::cout << "A root is: " << x << std::endl;

    std::cout << "Finding the root of \\sin(x/2) - 1 with initial values of 10 and 11" << std::endl;
    Secant_Method solver2(F1(), 10, 11);
    double x2 = solver2.solve();
    std::cout << "A root is: " << x2 << std::endl;
}

void solve2() {
    std::cout << "Finding the root of e^x - \\tanx with initial values of 1.0 and 1.4" << std::endl;
    Secant_Method solver(F2(), 1, 1.4);
    double x = solver.solve();
    std::cout << "A root is: " << x << std::endl;

    std::cout << "Finding the root of e^x - \\tanx with initial values of 2.3 and 2.4" << std::endl;
    Secant_Method solver2(F2(), 2.3, 2.4);
    double x2 = solver2.solve();
    std::cout << "A root is: " << x2 << std::endl;
}

void solve3() {
    std::cout << "Finding the root of x^3 - 12x^2 + 3x + 1 with initial values of 0 and -0.5" << std::endl;
    Secant_Method solver(F3(), 0, -0.5);
    double x = solver.solve();
    std::cout << "A root is: " << x << std::endl;

    std::cout << "Finding the root of x^3 - 12x^2 + 3x + 1 with initial values of 9 and 10" << std::endl;
    Secant_Method solver2(F3(), 9, 10);
    double x2 = solver2.solve();
    std::cout << "A root is: " << x2 << std::endl;
}

int main() {
    solve1();
    solve2();
    solve3();
    return 0;
}