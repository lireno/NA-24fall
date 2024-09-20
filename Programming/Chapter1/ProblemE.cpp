#include "EquationSolver.hpp"
#include "Function.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

const double Pi = acos(-1.);

class Volume : public Function {
  private:
    double L, V, r;

  public:
    Volume(double L, double V, double r) : L(L), V(V), r(r) {}
    double operator()(double h) const {
        return (0.5 * Pi * pow(r, 2) - pow(r, 2) * asin(h / r) - h * sqrt(pow(r, 2) - pow(h, 2))) * L - V;
        ;
    }
};

void solve() {
    double L = 10.0, r = 1.0, V = 12.4;
    Volume f(L, V, r);
    std::cout << "Finding the height with Bisection Method" << std::endl;
    Bisection_Method solver1(f, 0, 1);
    double h1 = solver1.solve();
    std::cout << "The height of the water in the tank is: " << std::fixed << std::setprecision(2) << h1 << std::endl;

    std::cout << "Finding the height with Newton Method" << std::endl;
    Newton_Method solver2(f, 0.3);
    double h2 = solver2.solve();
    std::cout << "The height of the water in the tank is: " << h2 << std::endl;

    std::cout << "Finding the height with Secant Method" << std::endl;
    Secant_Method solver3(f, 0.1, 0.2);
    double h3 = solver3.solve();
    std::cout << "The height of the water in the tank is: " << h3 << std::endl;
}

int main() {
    solve();
    return 0;
}