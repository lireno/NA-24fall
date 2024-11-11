// PPSpline.h
#ifndef PPSPLINE_H
#define PPSPLINE_H

#include "Function.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

class PPSpline : public Function {
  public:
    PPSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : nodes_(nodes), values_(values) {}

    PPSpline(const Function& f, const std::vector<double>& nodes)
        : nodes_(nodes) {
        for (double x : nodes) {
            values_.push_back(f(x));
        }
    }

    virtual ~PPSpline() = default;

    double evaluate(double x) const {
        if (nodes_.empty() || coefficients_.empty()) {
            throw std::runtime_error("Spline has not initialed.");
        }

        if (x < nodes_.front() || x > nodes_.back()) {
            throw std::out_of_range("Input x is out of the spline range.");
        }
        size_t interval = findInterval(x);

        const std::vector<double>& coeff = coefficients_[interval];

        double result = 0.0;
        double delta = x - nodes_[interval];
        for (auto it = coeff.rbegin(); it != coeff.rend(); ++it) {
            result = result * delta + *it;
        }

        return result;
    }

    double operator()(double x) const override {
        return evaluate(x);
    };

    void plotSpline(std::ostream& os = std::cout, size_t numSamplesPerInterval = 100) const {
        if (nodes_.empty() || coefficients_.empty()) {
            throw std::runtime_error("Spline not computed. Please call computeSpline() first.");
        }

        for (size_t i = 0; i < nodes_.size() - 1; ++i) {
            double xStart = nodes_[i];
            double xEnd = nodes_[i + 1];
            double step = (xEnd - xStart) / numSamplesPerInterval;

            for (size_t j = 0; j <= numSamplesPerInterval; ++j) {
                double x = xStart + j * step;
                double y = evaluate(x);
                os << x << " " << y << "\n";
            }
        }
    }

  protected:
    std::vector<double> nodes_;
    std::vector<double> values_;
    std::vector<std::vector<double>> coefficients_;

    virtual void computeSpline() = 0;

  private:
    // Find the interval that x belongs to
    size_t findInterval(double x) const {
        for (size_t i = 0; i < nodes_.size() - 1; ++i) {
            if (x >= nodes_[i] && x <= nodes_[i + 1]) {
                return i;
            }
        }
        throw std::runtime_error("Interval not found for input x.");
    }
};

class LinearPPSpline : public PPSpline {
  public:
    LinearPPSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : PPSpline(nodes, values) {
        computeSpline();
    }

    LinearPPSpline(const Function& f, const std::vector<double>& nodes)
        : PPSpline(f, nodes) {
        computeSpline();
    }

  protected:
    virtual void computeSpline() override {
        if (nodes_.size() != values_.size()) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }

        coefficients_.clear();
        for (size_t i = 0; i < nodes_.size() - 1; ++i) {
            double x1 = nodes_[i];
            double x2 = nodes_[i + 1];
            double y1 = values_[i];
            double y2 = values_[i + 1];

            double a = (y2 - y1) / (x2 - x1);
            double b = y1;

            coefficients_.push_back({b, a});
        }
    }
};

class CubicPPSpline : public PPSpline {
  public:
    CubicPPSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : PPSpline(nodes, values) {
        computePara();
    }

    CubicPPSpline(const Function& f, const std::vector<double>& nodes)
        : PPSpline(f, nodes) {
        computePara();
    }

  protected:
    std::vector<double> lambdas;
    std::vector<double> mus;
    std::vector<double> Ks; // Ks: divided differences

    // a: sub-diagonal, b: main diagonal, c: super-diagonal, d: right-hand side
    void thomasAlgorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x) {
        // cope with the case where b[0] is zero, which is needed in cyclicthomasAlgorithm
        if (std::abs(b[0]) < 1e-8) {
            std::vector<double> d_star(d.begin() + 1, d.end());
            double x2 = d[0] / c[0];
            d_star[0] -= x2 * b[1];
            d_star[1] -= x2 * a[1];
            std::vector<double> a_star(a.begin() + 1, a.end());
            a_star[0] = 0;
            std::vector<double> b_star(b.begin() + 1, b.end());
            b_star[0] = a[0];
            std::vector<double> c_star(c.begin() + 1, c.end());
            thomasAlgorithm(a_star, b_star, c_star, d_star, x);
            x.insert(x.begin() + 1, x2);
            x.pop_back();
            return;
        }
        int n = b.size();
        std::vector<double> c_prime(n - 1, 0.0);
        std::vector<double> d_prime(n, 0.0);

        // Forward elimination
        c_prime[0] = c[0] / b[0];
        d_prime[0] = d[0] / b[0];

        for (int i = 1; i < n - 1; ++i) {
            double m = b[i] - a[i - 1] * c_prime[i - 1];
            c_prime[i] = c[i] / m;
            d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) / m;
        }
        d_prime[n - 1] = (d[n - 1] - a[n - 2] * d_prime[n - 2]) / (b[n - 1] - a[n - 2] * c_prime[n - 2]);

        // Back substitution
        x[n - 1] = d_prime[n - 1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
        }
    }

    void cyclicthomasAlgorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x, const double an, const double cn) {
        int n = b.size();
        std::vector<double> u(n, 0.0);
        std::vector<double> v(n, 0.0);
        u[0] = b[0];
        v[0] = 1;
        u[n - 1] = cn;
        v[n - 1] = an / b[0];
        std::vector<double> b_star = b;
        b_star[n - 1] -= an * cn / b_star[0];
        b_star[0] = 0;
        std::vector<double> y(n, 0.0);
        std::vector<double> q(n, 0.0);
        thomasAlgorithm(a, b_star, c, d, y);
        thomasAlgorithm(a, b_star, c, u, q);
        double z = (v[0] * y[0] + v[n - 1] * y[n - 1]) / (1 + q[0] * v[0] + q[n - 1] * v[n - 1]);
        for (int i = 0; i < n; ++i) {
            x[i] = y[i] - z * q[i];
        }
    }

  private:
    void computePara() {
        lambdas.push_back(0);
        mus.push_back(0); // make the index consistent with the textbook
        for (size_t i = 1; i < nodes_.size() - 1; ++i) {
            double lambda = (nodes_[i + 1] - nodes_[i]) / (nodes_[i + 1] - nodes_[i - 1]);
            double mu = 1 - lambda;
            lambdas.push_back(lambda);
            mus.push_back(mu);
        }
        for (size_t i = 0; i < nodes_.size() - 1; ++i) {
            double dividedDifference = (values_[i + 1] - values_[i]) / (nodes_[i + 1] - nodes_[i]);
            Ks.push_back(dividedDifference);
        }
    }
};

class CompleteCubicPPSpline : public CubicPPSpline {
  public:
    CompleteCubicPPSpline(const std::vector<double>& nodes, const std::vector<double>& values, double da, double db)
        : CubicPPSpline(nodes, values) {
        derivative_a = da;
        derivative_b = db;
        computeSpline();
    }

    CompleteCubicPPSpline(const Function& f, const std::vector<double>& nodes)
        : CubicPPSpline(f, nodes) {
        computeDerivative(f);
        computeSpline();
    }

  protected:
    virtual void computeSpline() override {
        if (nodes_.size() != values_.size()) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }

        coefficients_.clear();

        size_t n = nodes_.size();
        std::vector<double> a(lambdas.begin() + 2, lambdas.begin() + n - 1); // Sub-diagonal
        std::vector<double> b(n - 2, 2.0);                                   // Main diagonal
        std::vector<double> c(mus.begin() + 1, mus.begin() + n - 2);         // Super-diagonal
        std::vector<double> ms(n - 2);
        std::vector<double> d(n - 2);
        for (size_t i = 1; i < n - 1; ++i) {
            d[i - 1] = 3 * (lambdas[i] * Ks[i - 1] + mus[i] * Ks[i]);
        }
        d[0] -= lambdas[1] * derivative_a;
        d[n - 3] -= mus[n - 2] * derivative_b;

        // Solve the system using Thomas algorithm
        thomasAlgorithm(a, b, c, d, ms);
        ms.insert(ms.begin(), derivative_a);
        ms.push_back(derivative_b);
        for (size_t i = 0; i < n - 1; ++i) {
            double c0 = values_[i];
            double c1 = ms[i];
            double h = nodes_[i + 1] - nodes_[i];
            double c2 = (3 * Ks[i] - 2 * ms[i] - ms[i + 1]) / h;
            double c3 = (ms[i + 1] + ms[i] - 2 * Ks[i]) / (h * h);
            coefficients_.push_back({c0, c1, c2, c3});
        }
    }

  private:
    double derivative_a;
    double derivative_b;
    void computeDerivative(const Function& f) {
        derivative_a = f.derivative(nodes_.front());
        derivative_b = f.derivative(nodes_.back());
    }
};

class NaturalCubicPPSpline : public CubicPPSpline {
  public:
    NaturalCubicPPSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : CubicPPSpline(nodes, values) {
        computeSpline();
    }

    NaturalCubicPPSpline(const Function& f, const std::vector<double>& nodes)
        : CubicPPSpline(f, nodes) {
        computeSpline();
    }

  protected:
    virtual void computeSpline() override {
        if (nodes_.size() != values_.size()) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }

        coefficients_.clear();

        size_t n = nodes_.size();
        std::vector<double> a(mus.begin() + 2, mus.begin() + n - 1);         // Sub-diagonal
        std::vector<double> b(n - 2, 2.0);                                   // Main diagonal
        std::vector<double> c(lambdas.begin() + 1, lambdas.begin() + n - 2); // Super-diagonal
        std::vector<double> Ms(n - 2);
        std::vector<double> d(n - 2);
        for (size_t i = 1; i < n - 1; ++i) {
            d[i - 1] = 3 * (Ks[i] - Ks[i - 1]) / (nodes_[i + 1] - nodes_[i - 1]);
        }
        thomasAlgorithm(a, b, c, d, Ms);
        Ms.insert(Ms.begin(), 0);
        Ms.push_back(0);
        for (size_t i = 0; i < n - 1; ++i) {
            double h = nodes_[i + 1] - nodes_[i];
            double c0 = values_[i];
            double c1 = Ks[i] - (Ms[i + 1] + 2 * Ms[i]) * h / 6;
            double c2 = Ms[i] / 2;
            double c3 = (Ms[i + 1] - Ms[i]) / (6 * h);
            coefficients_.push_back({c0, c1, c2, c3});
        }
    }
};

class PeriodicCubicPPSpline : public CubicPPSpline {
  public:
    PeriodicCubicPPSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : CubicPPSpline(nodes, values) {
        if (std::abs(values.front() - values.back()) > 1e-7) {
            throw std::runtime_error("The first and last values must be the same for a periodic spline.");
        }
        computeSpline();
    }

    PeriodicCubicPPSpline(const Function& f, const std::vector<double>& nodes)
        : CubicPPSpline(f, nodes) {
        computeSpline();
    }

  protected:
    virtual void computeSpline() override {
        if (nodes_.size() != values_.size()) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }

        coefficients_.clear();

        size_t n = nodes_.size();
        std::vector<double> a(lambdas.begin() + 2, lambdas.begin() + n - 1); // Sub-diagonal
        std::vector<double> b(n - 2, 2.0);                                   // Main diagonal
        std::vector<double> c(mus.begin() + 1, mus.begin() + n - 2);         // Super-diagonal
        std::vector<double> ms(n - 2);
        std::vector<double> d(n - 2);
        for (size_t i = 1; i < n - 1; ++i) {
            d[i - 1] = 3 * (lambdas[i] * Ks[i - 1] + mus[i] * Ks[i]);
        }
        double h1 = nodes_[1] - nodes_[0];
        double h2 = nodes_[n - 1] - nodes_[n - 2];
        double temp = 1.5 * (Ks[0] * h2 + Ks[n - 2] * h1) / (h1 + h2);
        d[0] -= lambdas[1] * temp;
        d[n - 3] -= mus[n - 2] * temp;
        b[0] -= lambdas[1] * h2 / (h1 + h2);
        b[n - 3] -= mus[n - 2] * h1 / (h1 + h2);
        double an = -lambdas[1] * 0.5 * h1 / (h1 + h2);
        double cn = -mus[n - 2] * 0.5 * h2 / (h1 + h2);
        cyclicthomasAlgorithm(a, b, c, d, ms, an, cn);
        double m1 = 2 * h1 / (h1 + h2) * ms[0];
        ms.push_back(m1);
        ms.insert(ms.begin(), m1);
        for (size_t i = 0; i < n - 1; ++i) {
            double c0 = values_[i];
            double c1 = ms[i];
            double h = nodes_[i + 1] - nodes_[i];
            double c2 = (3 * Ks[i] - 2 * ms[i] - ms[i + 1]) / h;
            double c3 = (ms[i + 1] + ms[i] - 2 * Ks[i]) / (h * h);
            coefficients_.push_back({c0, c1, c2, c3});
        }
    }
};
#endif // PPSPLINE_H
