// PPSpline.h
// This file defines classes for constructing and evaluating piecewise polynomial splines,
// including linear and cubic variations. The splines are used to interpolate a set of
// nodes and values or a given function. The main classes are:

// - PPSpline: Abstract base class for general piecewise polynomial splines.
// - LinearPPSpline: Implements a linear interpolation spline.
// - CubicPPSpline: Base class for cubic splines, with subclasses for:
//      - CompleteCubicPPSpline: Specifies derivatives at boundaries.
//      - NaturalCubicPPSpline: Sets the second derivatives at boundaries to zero.
//      - PeriodicCubicPPSpline: Ensures periodicity in the spline, with equal values at the boundaries.

// Key methods include:
// - evaluate(double x): Computes the interpolated value at a given point.
// - computeSpline(): Virtual method implemented by derived classes to compute spline coefficients.
// - plot(): Generates a plot of the spline for visualization.

#ifndef PPSPLINE_H
#define PPSPLINE_H

#include "Function.h"
#include "Thomas.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

class PPSpline : public Function {
  public:
    PPSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : nodes_(nodes), values_(values) { n = nodes_.size(); }

    PPSpline(const Function& f, const std::vector<double>& nodes)
        : nodes_(nodes) {
        for (double x : nodes) {
            values_.push_back(f(x));
        }
        n = nodes_.size();
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

    void plot(std::ostream& os = std::cout, size_t numSample = 300) const {
        Function::plot(os, nodes_.front(), nodes_.back(), numSample);
    }

  protected:
    std::vector<double> nodes_;
    std::vector<double> values_;
    std::vector<std::vector<double>> coefficients_;
    int n; // size of nodes_

    virtual void computeSpline() = 0;

  private:
    // Find the interval that x belongs to
    size_t findInterval(double x) const {
        if (x < nodes_.front() || x > nodes_.back()) {
            throw std::out_of_range("Input x is out of the spline range.");
        }

        size_t low = 0;
        size_t high = n - 1;

        while (low < high) {
            size_t mid = low + (high - low) / 2;
            if (x < nodes_[mid]) {
                high = mid;
            } else if (x > nodes_[mid + 1]) {
                low = mid + 1;
            } else {
                return mid;
            }
        }

        return low;
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
        for (size_t i = 0; i < n - 1; ++i) {
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

  private:
    void computePara() {
        lambdas.push_back(0);
        mus.push_back(0); // make the index consistent with the textbook
        for (size_t i = 1; i < n - 1; ++i) {
            double lambda = (nodes_[i + 1] - nodes_[i]) / (nodes_[i + 1] - nodes_[i - 1]);
            double mu = 1 - lambda;
            lambdas.push_back(lambda);
            mus.push_back(mu);
        }
        for (size_t i = 0; i < n - 1; ++i) {
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

        std::vector<double> a(mus.begin() + 2, mus.begin() + n - 1);         // Sub-diagonal
        std::vector<double> b(n - 2, 2.0);                                   // Main diagonal
        std::vector<double> c(lambdas.begin() + 1, lambdas.begin() + n - 2); // Super-diagonal
        std::vector<double> Ms(n - 2);
        std::vector<double> d(n - 2);
        for (size_t i = 1; i < n - 1; ++i) {
            d[i - 1] = 6 * (Ks[i] - Ks[i - 1]) / (nodes_[i + 1] - nodes_[i - 1]);
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

        double h1 = nodes_[1] - nodes_[0];
        double h2 = nodes_[n - 1] - nodes_[n - 2];
        lambdas[0] = h1 / (h1 + h2);
        mus[0] = h2 / (h1 + h2);

        std::vector<double> d(n - 1);
        d[0] = 3 * (lambdas[0] * Ks[n - 2] + mus[0] * Ks[0]);
        for (size_t i = 1; i < n - 1; ++i) {
            d[i] = 3 * (lambdas[i] * Ks[i - 1] + mus[i] * Ks[i]);
        }

        std::vector<double> b(n - 1, 2.0); // Main diagonal
        std::vector<double> ms(n - 1);

        cyclicthomasAlgorithm(lambdas, b, mus, d, ms);
        ms.push_back(ms.front());

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
