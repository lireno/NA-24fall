// BSpline.h
// This file defines classes for B-spline interpolation and evaluation, including support for linear, cubic, and advanced variations.
// These classes can interpolate data points or a given function, using B-spline basis functions and coefficients. Key classes include:

// - Bbase: Represents a single B-spline basis function for a given set of nodes.
//      - Methods include evaluation, first derivative, and second derivative computation of the basis function.

// - BSpline: Abstract base class for general B-spline interpolation.
//      - Supports initialization with nodes and values or a Function object.
//      - Key methods include evaluate() for computing spline values and plot() for visualization.

// Derived classes:
// - LinearBSpline: Implements linear B-spline interpolation.
// - QuarticBSpline: Implements quartic B-spline interpolation with support for extended nodes and basis functions.
// - CubicBSpline: Implements cubic B-spline interpolation with support for extended nodes and basis functions.
//      - Variations include:
//          - CompleteCubicBSpline: Specifies boundary derivatives for interpolation.
//          - NaturalCubicBSpline: Ensures second derivatives are zero at the boundaries.
//          - PeriodicCubicBSpline: Ensures periodic boundary conditions for seamless looping.

#ifndef BSPLINE_H
#define BSPLINE_H

#include "Function.h"
#include "Thomas.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

class Bbase : public Function {
  public:
    Bbase(const std::vector<double>& nodes)
        : nodes_(nodes) { degree_ = nodes_.size() - 2; }

    double evaluate(double x) const {
        if (nodes_.empty()) {
            throw std::runtime_error("Spline has not initialed.");
        }
        if (x <= nodes_.front() || x >= nodes_.back()) {
            return 0;
        }
        if (degree_ == 0) {
            return 1;
        }
        if (degree_ == 1) {
            if (x > nodes_[1]) {
                return (nodes_[2] - x) / (nodes_[2] - nodes_[1]);
            } else {
                return (x - nodes_[0]) / (nodes_[1] - nodes_[0]);
            }
        } else {
            std::vector<double> nodes1(nodes_.begin(), nodes_.end() - 1);
            Bbase b1(nodes1);

            std::vector<double> nodes2(nodes_.begin() + 1, nodes_.end());
            Bbase b2(nodes2);

            double result = 0;
            result += (x - nodes_.front()) / (nodes_[degree_] - nodes_.front()) * b1.evaluate(x);
            result += (nodes_.back() - x) / (nodes_.back() - nodes_[1]) * b2.evaluate(x);

            return result;
        }
    }
    double operator()(double x) const override {
        return evaluate(x);
    }

    double derivative(double x) const override {
        if (nodes_.empty()) {
            throw std::runtime_error("Spline has not initialed.");
        }
        if (x <= nodes_.front() || x >= nodes_.back()) {
            return 0;
        }
        if (degree_ == 0) {
            return 0;
        } else {
            std::vector<double> nodes1(nodes_.begin(), nodes_.end() - 1);
            Bbase b1(nodes1);

            std::vector<double> nodes2(nodes_.begin() + 1, nodes_.end());
            Bbase b2(nodes2);

            double result = 0;
            result += b1.evaluate(x) / (nodes_[degree_] - nodes_.front());
            result -= b2.evaluate(x) / (nodes_.back() - nodes_[1]);
            result *= degree_;

            return result;
        }
    }

    double secondDerivative(double x) const {
        if (nodes_.empty()) {
            throw std::runtime_error("Spline has not initialed.");
        }
        if (x <= nodes_.front() || x >= nodes_.back()) {
            return 0;
        }
        if (degree_ == 0) {
            return 0;
        } else {
            std::vector<double> nodes1(nodes_.begin(), nodes_.end() - 1);
            Bbase b1(nodes1);

            std::vector<double> nodes2(nodes_.begin() + 1, nodes_.end());
            Bbase b2(nodes2);

            double result = 0;
            result += b1.derivative(x) / (nodes_[degree_] - nodes_.front());
            result -= b2.derivative(x) / (nodes_.back() - nodes_[1]);
            result *= degree_;

            return result;
        }
    }

  private:
    std::vector<double> nodes_;
    int degree_;
};

class BSpline : public Function {
  public:
    BSpline() = default;

    BSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : nodes_(nodes), values_(values) { n = nodes.size(); }

    BSpline(const Function& f, const std::vector<double>& nodes)
        : nodes_(nodes) {
        for (double x : nodes) {
            values_.push_back(f(x));
        }
        n = nodes.size();
    }

    BSpline(const std::vector<double>& nodes, const std::vector<double>& coefficients, int degree)
        : nodes_(nodes), coefficients_(coefficients), degree_(degree) {
        n = nodes.size();
        int N = coefficients.size();
        if (N != n + degree - 1) {
            throw std::runtime_error("Number of coefficients and nodes are not fit with degree.");
        }
        generate_extended_points();
        generate_basis();
    }

    double evaluate(double x) const {
        if (nodes_.empty() || coefficients_.empty()) {
            throw std::runtime_error("Spline has not initialed.");
        }

        double result = 0.0;

        for (size_t i = 0; i < coefficients_.size(); ++i) {
            result = result + basis_[i](x) * coefficients_[i];
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
    std::vector<double> extended_nodes_;
    std::vector<double> values_;
    std::vector<double> coefficients_;
    std::vector<Bbase> basis_;
    int degree_;
    int n; // size of nodes_

    void generate_extended_points() {
        int degree = degree_;

        extended_nodes_.clear();
        double leftExtension = nodes_[1] - nodes_[0];
        double rightExtension = nodes_[n - 1] - nodes_[n - 2];

        for (int i = 0; i < degree; ++i) {
            extended_nodes_.push_back(nodes_[0] - (degree - i) * leftExtension);
        }
        for (double node : nodes_) {
            extended_nodes_.push_back(node);
        }
        for (int i = 0; i < degree; ++i) {
            extended_nodes_.push_back(nodes_[n - 1] + (i + 1) * rightExtension);
        }
    }

    void generate_periodic_extended_points() {
        extended_nodes_.clear();
        int degree = degree_;

        for (int i = 0; i < degree; ++i) {
            extended_nodes_.push_back(nodes_.front() - (nodes_.back() - nodes_[n - degree + i - 1]));
        }
        for (double node : nodes_) {
            extended_nodes_.push_back(node);
        }
        for (int i = 0; i < degree; ++i) {
            extended_nodes_.push_back(nodes_.back() + (nodes_[i + 1] - nodes_.front()));
        }
    }

    void generate_basis() {
        basis_.clear();
        int numBasis = n + degree_ - 1;
        for (int i = 0; i < numBasis; ++i) {
            std::vector<double> basisNodes(extended_nodes_.begin() + i, extended_nodes_.begin() + i + degree_ + 2);
            Bbase basisFunction(basisNodes);
            basis_.push_back(basisFunction);
        }
    }
};

class LinearBSpline : public BSpline {
  public:
    LinearBSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : BSpline(nodes, values) {
        degree_ = 1;
        generate_extended_points();
        generate_basis();
        computeSpline();
    }

    LinearBSpline(const Function& f, const std::vector<double>& nodes)
        : BSpline(f, nodes) {
        degree_ = 1;
        generate_extended_points();
        generate_basis();
        computeSpline();
    }

  protected:
    void computeSpline() {
        coefficients_ = values_;
    }
};

class quarticBSpline : public BSpline {
  public:
    quarticBSpline(const Function& f, const int l, const int r)
        : f(f), l(l), r(r) {
        degree_ = 2;
        n = r - l + 1;
        for (int i = l; i <= r; ++i) {
            nodes_.push_back(i);
        }
        generate_extended_points();
        generate_basis();
        computeSpline();
    };

    void generate_extended_points() {
        extended_nodes_.clear();
        for (int i = l - 2; i <= r + 2; ++i) {
            extended_nodes_.push_back(i);
        }
    }

    void computeSpline() {
        coefficients_.clear();
        std::vector<double> a(n - 2, 1.0);
        std::vector<double> b(n - 1, 6.0);
        std::vector<double> c(n - 2, 1.0);
        std::vector<double> d(n - 1, 0.0);
        b.front() = 5.0;
        b.back() = 5.0;
        for (int i = 0; i < n - 1; ++i) {
            d[i] = 8 * f(l + i + 1.0 / 2.0);
        }
        d[0] -= 2 * f(l);
        d[n - 2] -= 2 * f(r);
        coefficients_.resize(n - 1);
        thomasAlgorithm(a, b, c, d, coefficients_);
        double a0 = 2 * f(l) - coefficients_.front();
        double an = 2 * f(r) - coefficients_.back();
        coefficients_.insert(coefficients_.begin(), a0);
        coefficients_.push_back(an);
    }

  private:
    int l;
    int r;
    const Function& f;
};

class CubicBSpline : public BSpline {
  public:
    CubicBSpline() = default;
    CubicBSpline(const std::vector<double>& nodes, const std::vector<double>& values, bool periodic = false)
        : BSpline(nodes, values) {
        degree_ = 3;
        if (periodic) {
            generate_periodic_extended_points();
        } else {
            generate_extended_points();
        }
        generate_basis();
        computeVal();
    }

    CubicBSpline(const Function& f, const std::vector<double>& nodes, bool periodic = false)
        : BSpline(f, nodes) {
        degree_ = 3;
        if (periodic) {
            generate_periodic_extended_points();
        } else {
            generate_extended_points();
        }
        generate_basis();
        computeVal();
    }

  protected:
    std::vector<double> a; // store the values of nodes for each basis
    std::vector<double> b;
    std::vector<double> c;

    void solveUniqueTridialog(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x, const std::vector<double>& s, const std::vector<double>& t) {
        std::vector<double> a_star = a;
        std::vector<double> b_star = b;
        std::vector<double> c_star = c;
        std::vector<double> d_star = d;
        double m1 = a_star[0] / s[0];
        b_star[0] = b_star[0] - m1 * s[1];
        c_star[0] = c_star[0] - m1 * s[2];
        d_star[0] = d_star[0] - m1 * s[3];
        double m2 = c[n - 1] / t[2];
        a_star[n - 1] = a_star[n - 1] - m2 * t[0];
        b_star[n - 1] = b_star[n - 1] - m2 * t[1];
        d_star[n - 1] = d_star[n - 1] - m2 * t[3];
        a_star.erase(a_star.begin());
        c_star.pop_back();
        x.resize(n);
        thomasAlgorithm(a_star, b_star, c_star, d_star, x);
        double a0 = (s[3] - s[1] * x[0] - s[2] * x[1]) / s[0];
        double an = (t[3] - t[0] * x[n - 2] - t[1] * x[n - 1]) / t[2];
        x.insert(x.begin(), a0);
        x.push_back(an);
    };

  private:
    void computeVal() {
        int n = nodes_.size();
        a.resize(n);
        b.resize(n);
        c.resize(n);
        for (int i = 0; i < n; ++i) {
            a[i] = basis_[i](nodes_[i]);
            b[i] = basis_[i + 1](nodes_[i]);
            c[i] = basis_[i + 2](nodes_[i]);
        }
    }
};

class CompleteCubicBSpline : public CubicBSpline {
  public:
    CompleteCubicBSpline(const std::vector<double>& nodes, const std::vector<double>& values, double da, double db)
        : CubicBSpline(nodes, values) {
        derivative_a = da;
        derivative_b = db;
        computeSpline();
    }

    CompleteCubicBSpline(const Function& f, const std::vector<double>& nodes)
        : CubicBSpline(f, nodes) {
        computeDerivative(f);
        computeSpline();
    }

  protected:
    void computeSpline() {
        if (values_.size() != n) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }
        coefficients_.clear();
        std::vector<double> s(4);
        std::vector<double> t(4);
        double x0 = nodes_[0];
        double xn = nodes_[n - 1];
        for (int i = 0; i < 3; ++i) {
            s[i] = basis_[i].derivative(x0);
            t[i] = basis_[n - 1 + i].derivative(xn);
        }
        s[3] = derivative_a;
        t[3] = derivative_b;
        solveUniqueTridialog(a, b, c, values_, coefficients_, s, t);
    };

  private:
    double derivative_a;
    double derivative_b;
    void computeDerivative(const Function& f) {
        derivative_a = f.derivative(nodes_.front());
        derivative_b = f.derivative(nodes_.back());
    }
};

class NaturalCubicBSpline : public CubicBSpline {
  public:
    NaturalCubicBSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : CubicBSpline(nodes, values) {
        computeSpline();
    }

    NaturalCubicBSpline(const Function& f, const std::vector<double>& nodes)
        : CubicBSpline(f, nodes) {
        computeSpline();
    }

  protected:
    void computeSpline() {
        if (values_.size() != n) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }
        coefficients_.clear();
        std::vector<double> s(4);
        std::vector<double> t(4);
        double x0 = nodes_[0];
        double xn = nodes_[n - 1];
        for (int i = 0; i < 3; ++i) {
            s[i] = basis_[i].secondDerivative(x0);
            t[i] = basis_[n - 1 + i].secondDerivative(xn);
        }
        s[3] = 0;
        t[3] = 0;
        solveUniqueTridialog(a, b, c, values_, coefficients_, s, t);
    };
};

class PeriodicCubicBSpline : public CubicBSpline {
  public:
    PeriodicCubicBSpline() = default;

    PeriodicCubicBSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : CubicBSpline(nodes, values, true) {
        if (std::abs(values.front() - values.back()) > 1e-7) {
            throw std::runtime_error("The first and last values must be the same for a periodic spline.");
        }
        computeSpline();
    }

    PeriodicCubicBSpline(const Function& f, const std::vector<double>& nodes)
        : CubicBSpline(f, nodes, true) {
        computeSpline();
    }

  protected:
    void computeSpline() {
        if (values_.size() != n) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }

        coefficients_.clear();
        std::vector<double> a_star(a);
        std::vector<double> b_star(b);
        std::vector<double> c_star(c);
        std::vector<double> d_star(values_);

        a_star.pop_back();
        b_star.pop_back();
        c_star.pop_back();
        d_star.pop_back();
        coefficients_.resize(n - 1);

        cyclicthomasAlgorithm(a_star, b_star, c_star, d_star, coefficients_);
        coefficients_.insert(coefficients_.begin(), coefficients_.back());
        coefficients_.push_back(coefficients_[1]);
        coefficients_.push_back(coefficients_[2]);
    };
};

#endif // BSPLINE_H