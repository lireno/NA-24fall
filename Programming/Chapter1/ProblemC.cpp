#include "EquationSolver.hpp"
#include "Function.hpp"
#include <cmath>
#include <iostream>

class F : public Function {
  public:
    double operator()(double x) const {
        return tan(x) - x;
    }
};

void solve1() {
    std::cout << "Finding the root of \\tan x - x near 4.5" << std::endl;
    Newton_Method solver(F(), 4.5);
    double x = solver.solve();
    std::cout << "The root near 4.5 is: " << x << std::endl;
}

void solve2() {
    std::cout << "Finding the root of \\tan x - x near 4.7" << std::endl;
    Newton_Method solver(F(), 4.7);
    double x = solver.solve();
    std::cout << "The root near 4.7 is: " << x << std::endl;
}

int main() {
    solve1();
    solve2();
    return 0;
}