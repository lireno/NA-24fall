// PPformSpline.h
#ifndef PPFORMSPLINE_H
#define PPFORMSPLINE_H

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

class Function {
  public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const {
        double h = 1e-6;
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h);
    }
};

class PPformSpline {
  public:
    PPformSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : nodes_(nodes), values_(values) {}

    PPformSpline(const Function& f, const std::vector<double>& nodes)
        : nodes_(nodes) {
        for (double x : nodes) {
            values_.push_back(f(x));
        }
    }

    virtual ~PPformSpline() = default;

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
    size_t findInterval(double x) const {
        for (size_t i = 0; i < nodes_.size() - 1; ++i) {
            if (x >= nodes_[i] && x <= nodes_[i + 1]) {
                return i;
            }
        }
        throw std::runtime_error("Interval not found for input x.");
    }
};

class LinearPPformSpline : public PPformSpline {
  public:
    LinearPPformSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : PPformSpline(nodes, values) {
        computeSpline();
    }

    LinearPPformSpline(const Function& f, const std::vector<double>& nodes)
        : PPformSpline(f, nodes) {
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

class CubicPPformSpline : public PPformSpline {
  public:
    CubicPPformSpline(const std::vector<double>& nodes, const std::vector<double>& values)
        : PPformSpline(nodes, values) {
        computePara();
    }

    CubicPPformSpline(const Function& f, const std::vector<double>& nodes)
        : PPformSpline(f, nodes) {
        computePara();
    }

  protected:
    std::vector<double> lambdas;
    std::vector<double> mus;
    std::vector<double> Ks; // Ks: divided differences

    // a: sub-diagonal, b: main diagonal, c: super-diagonal, d: right-hand side
    void thomasAlgorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x) {
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

class CompleteCubicPPformSPline : public CubicPPformSpline {
  public:
    CompleteCubicPPformSPline(const std::vector<double>& nodes, const std::vector<double>& values, double da, double db)
        : CubicPPformSpline(nodes, values) {
        derivative_a = da;
        derivative_b = db;
        computeSpline();
    }

    CompleteCubicPPformSPline(const Function& f, const std::vector<double>& nodes)
        : CubicPPformSpline(f, nodes) {
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
        for (const auto& coeff : coefficients_) {
            std::cout << "系数: ";
            for (double c : coeff) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
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

class NaturalCubicPPformSPline : public CubicPPformSpline {
  public:
    NaturalCubicPPformSPline(const std::vector<double>& nodes, const std::vector<double>& values, double derivative_a, double derivative_b)
        : CubicPPformSpline(nodes, values) {
        computeSpline();
    }

    NaturalCubicPPformSPline(const Function& f, const std::vector<double>& nodes)
        : CubicPPformSpline(f, nodes) {

        computeSpline();
    }

  protected:
    virtual void computeSpline() override {
        if (nodes_.size() != values_.size()) {
            throw std::runtime_error("Number of nodes and values must be the same.");
        }

        coefficients_.clear();
    }
};

class PeriodicCubicPPformSPline : public CubicPPformSpline {
  public:
    PeriodicCubicPPformSPline(const std::vector<double>& nodes, const std::vector<double>& values)
        : CubicPPformSpline(nodes, values) {
        if (std::abs(values.front() - values.back()) > 1e-7) {
            throw std::runtime_error("The first and last values must be the same for a periodic spline.");
        }
        computeSpline();
    }

    PeriodicCubicPPformSPline(const Function& f, const std::vector<double>& nodes)
        : CubicPPformSpline(f, nodes) {
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
        double an = -lambdas[1] * h1 / (h1 + h2);
        double cn = -mus[n - 2] * h2 / (h1 + h2);
        cyclicthomasAlgorithm(a, b, c, d, ms, an, cn);
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
#endif // PPFORMSPLINE_H
