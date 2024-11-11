// PPSpline.h
#ifndef BSPLINE_H
#define BSPLINE_H

#include "Function.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

class BSpline : public Function {
  public:
    BSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : nodes_(nodes), values_(values) {}

    BSpline(const Function& f, const std::vector<double>& nodes)
        : nodes_(nodes) {
        for (double x : nodes) {
            values_.push_back(f(x));
        }
    }

    virtual ~BSpline() = default;

    double evaluate(double x) const {
        if (nodes_.empty() || coefficients_.empty()) {
            throw std::runtime_error("Spline has not initialed.");
        }

        if (x < nodes_.front() || x > nodes_.back()) {
            throw std::out_of_range("Input x is out of the spline range.");
        }

        double result = 0.0;

        for (size_t i = 0; i < coefficients_.size() - 1; ++i) {
            result = result + basis_[i](x) * coefficients_[i];
        }

        return result;
    }

    double operator()(double x) const override {
        return evaluate(x);
    };

    void plotSpline(std::ostream& os = std::cout, size_t numSamplesPerInterval = 100) const {
        if (nodes_.empty()) {
            throw std::runtime_error("Spline not computed.");
        }

        for (size_t i = 0; i < nodes_.size() - 1; ++i) {
            double xStart = nodes_[i];
            double xEnd = nodes_[i + 1];
            double step = (xEnd - xStart) / numSamplesPerInterval;

            for (size_t j = 0; j <= numSamplesPerInterval; ++j) {
                double x = xStart + j * step;
                os << x << " " << evaluate(x) << std::endl;
            }
        }
    }

  protected:
    std::vector<double> nodes_;
    std::vector<double> values_;
    std::vector<double> coefficients_;
    std::vector<Bbase> basis_;
    virtual void computeSpline() = 0;
};

class Bbase : public Function {
  public:
    Bbase(const std::vector<double>& nodes, const std::vector<double>& values)
        : nodes_(nodes), values_(values) { degree_ = nodes_.size() - 2; }

    Bbase(const Function& f, const std::vector<double>& nodes)
        : nodes_(nodes) {
        for (double x : nodes) {
            values_.push_back(f(x));
        }
        degree_ = nodes_.size() - 2;
    }

    double evaluate(double x) const {
        if (nodes_.empty()) {
            throw std::runtime_error("Spline has not initialed.");
        }
        if (x < nodes_.front() || x > nodes_.back()) {
            return 0;
        }
        if (degree_ == 0) {
            return 1;
        }
        if (degree_ == 1) {
            if (x > nodes_[1]) {
                return (nodes_[1] - x) / (nodes_[1] - nodes_[0]);
            } else {
                return (x - nodes_[0]) / (nodes_[1] - nodes_[0]);
            }
        } else {
            std::vector<double> nodes1(nodes_.begin(), nodes_.end() - 1);
            std::vector<double> values1(values_.begin(), values_.end() - 1);
            Bbase b1(nodes1, values1);

            std::vector<double> nodes2(nodes_.begin() + 1, nodes_.end());
            std::vector<double> values2(values_.begin() + 1, values_.end());
            Bbase b2(nodes2, values2);

            double result = 0;
            result += (x - nodes_.front()) / (nodes_[degree_] - nodes_.front()) * b1.evaluate(x);
            result += (nodes_.back() - x) / (nodes_.back() - nodes_[1]) * b2.evaluate(x);

            return result;
        }
    }

  public:
    double operator()(double x) const override {
        return evaluate(x);
    }

  private:
    std::vector<double> nodes_;
    std::vector<double> values_;
    int degree_;
};

#endif // BSPLINE_H