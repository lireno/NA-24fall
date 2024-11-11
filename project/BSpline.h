// PPSpline.h
#ifndef BSPLINE_H
#define BSPLINE_H

#include "Function.h"
#include "Thomas.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
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
                return (nodes_[2] - x) / (nodes_[1] - nodes_[0]);
            } else {
                return (x - nodes_[0]) / (nodes_[2] - nodes_[1]);
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

  private:
    std::vector<double> nodes_;
    int degree_;
};

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

        double result = 0.0;

        for (size_t i = 0; i < coefficients_.size() + degree_ - 1; ++i) {
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
    std::vector<double> extended_nodes_;
    std::vector<double> values_;
    std::vector<double> coefficients_;
    std::vector<Bbase> basis_;
    int degree_;
    virtual void computeSpline() = 0;
    void generatebasis() {
        int n = nodes_.size();
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
        basis_.clear();
        int numBasis = n + degree - 1;
        for (int i = 0; i < numBasis; ++i) {
            std::vector<double> basisNodes(extended_nodes_.begin() + i, extended_nodes_.begin() + i + degree + 2);
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
        generatebasis();
        computeSpline();
    }

    LinearBSpline(const Function& f, const std::vector<double>& nodes)
        : BSpline(f, nodes) {
        degree_ = 1;
        generatebasis();
        computeSpline();
    }

  protected:
    virtual void computeSpline() override {
        coefficients_ = values_;
    }
};

class CubicBSpline : public BSpline {
  public:
    CubicBSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : BSpline(nodes, values) {
        degree_ = 3;
        generatebasis();
        computeSpline();
    }

    CubicBSpline(const Function& f, const std::vector<double>& nodes)
        : BSpline(f, nodes) {
        degree_ = 3;
        generatebasis();
        computeSpline();
    }

  protected:
    virtual void computeSpline() override {
        if (nodes_.size() != values_.size()) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }
    }
};

#endif // BSPLINE_H