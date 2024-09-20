#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER
#include "Function.hpp"
#include <cmath>
#include <stdexcept>

class EquationSolver {
  protected:
    const Function& F;

  public:
    EquationSolver(const Function& F) : F(F) {}
    virtual double solve() = 0;
};

class Bisection_Method : public EquationSolver {
  private:
    double a, b;
    double eps, delta;
    int Maxiter;

  public:
    Bisection_Method(const Function& F, double a, double b,
                     double eps = 1e-7, double delta = 1e-6, int Maxiter = 50)
        : EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}

    virtual double solve() {
        double fa = F(a), fb = F(b);

        // Check if the function changes sign over the interval
        if (fa * fb > 0) {
            throw std::invalid_argument("The function must change sign over the interval.");
        }

        for (int iter = 0; iter < Maxiter; ++iter) {
            double c = (a + b) / 2;
            double fc = F(c);

            // Check if the function is continuous in the interval
            if (F.is_Continuous(c) == false) {
                throw std::runtime_error("The function appears to be discontinuous in the interval.");
            }

            // Check if the precision is met
            if (fabs(fc) < delta || (b - a) / 2 < eps) {
                return c;
            }

            // adjust the interval
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }
        }

        return (a + b) / 2; // return the best guess after max iterations
    }
};

class Newton_Method : public EquationSolver {
  private:
    double x0;
    double eps;
    int Maxiter;

  public:
    Newton_Method(const Function& F, double x0,
                  double eps = 1e-7, double Maxiter = 8) : EquationSolver(F), x0(x0), Maxiter(Maxiter), eps(eps) {}

    virtual double solve() {

        double x = x0;

        for (int iter = 0; iter < Maxiter; ++iter) {
            double fx = F(x);
            double dfx = F.derivative(x);

            if (fabs(fx) < eps) {
                return x;
            }

            if (fabs(dfx) < 1e-12) { // avoid division by zero
                throw std::runtime_error("Derivative too small.");
            }

            x = x - fx / dfx;
        }

        return x; // return the best guess after max iterations
    }
};

class Secant_Method : public EquationSolver {
  private:
    double x0, x1;
    double eps;
    int Maxiter;

  public:
    Secant_Method(const Function& F, double x0, double x1,
                  double eps = 1e-7, int Maxiter = 50) : EquationSolver(F), x0(x0), x1(x1), Maxiter(Maxiter), eps(eps) {}

    virtual double solve() {
        double x_prev = x0, x_curr = x1;

        for (int iter = 0; iter < Maxiter; ++iter) {
            double f_prev = F(x_prev);
            double f_curr = F(x_curr);

            if (fabs(f_curr) < eps) {
                return x_curr;
            }

            if (fabs(f_curr - f_prev) < 1e-12) { // avoid division by zero
                throw std::runtime_error("Function values too close.");
            }

            double x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);

            x_prev = x_curr;
            x_curr = x_next;
        }

        return x_curr; // return the best guess after max iterations
    }
};

#endif