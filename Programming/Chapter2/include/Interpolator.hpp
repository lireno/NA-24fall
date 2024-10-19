#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP

#include "Function.hpp"
#include <iostream>
#include <vector>

class Interpolator {
  public:
    virtual double interpolate(double x) const = 0;
    virtual ~Interpolator() = default;
};

class point {
  public:
    point() : x(0), y(0) {}
    point(double x, double y) : x(x), y(y) {}

    point operator+(const point& other) const {
        return point(x + other.x, y + other.y);
    }

    point operator-(const point& other) const {
        return point(x - other.x, y - other.y);
    }

    point operator*(double scalar) const {
        return point(x * scalar, y * scalar);
    }

    point operator/(double scalar) const {
        return point(x / scalar, y / scalar);
    }

    void normalize() {
        double length = std::sqrt(x * x + y * y);
        x /= length;
        y /= length;
    }

    void print() const {
        std::cout << "(" << x << ", " << y << ")";
    }

    double getX() const { return x; }
    double getY() const { return y; }

  private:
    double x;
    double y;
};

class PointInterpolator {
  public:
    virtual point Interpolate(double t) const = 0;
    virtual ~PointInterpolator() = default;
};

class NewtonInterpolator : public Interpolator {
  public:
    NewtonInterpolator(const std::vector<double>& x_values, const std::vector<double>& f_values)
        : x_values(x_values), f_values(f_values), divided_differences(x_values.size(), std::vector<double>(x_values.size(), 0.0)) {
        if (x_values.empty() || f_values.empty() || x_values.size() != f_values.size()) {
            throw std::invalid_argument("x_values and f_values must be non-empty and of the same size");
        }
        computeDividedDifferences();
    }

    NewtonInterpolator(const Function& func, const std::vector<double>& x_values)
        : x_values(x_values), divided_differences(x_values.size(), std::vector<double>(x_values.size(), 0.0)) {
        if (x_values.empty()) {
            throw std::invalid_argument("x_values cannot be empty");
        }
        f_values.resize(x_values.size());
        for (int i = 0; i < x_values.size(); ++i) {
            f_values[i] = func(x_values[i]);
        }
        computeDividedDifferences();
    }

    double interpolate(double x_value) const override {
        int n = x_values.size();
        double result = divided_differences[n - 1][n - 1];

        for (int i = n - 1; i >= 0; --i) {
            result = divided_differences[i][i] + (x_value - x_values[i]) * result;
        }

        return result;
    }

    void printInternalData() const {
        std::cout << "x-values: ";
        for (const auto& val : x_values) {
            std::cout << val << " ";
        }
        std::cout << "\nDivided differences: ";

        for (int i = 0; i < divided_differences.size(); ++i) {
            std::cout << divided_differences[i][i] << " ";
        }
        std::cout << "\n";

        // the following code is for debugging purposes
        // std::cout << "\nDivided differences table: \n";
        // for (int i = 0; i < divided_differences.size(); ++i) {
        //     for (int j = 0; j < i + 1; ++j) {
        //         std::cout << divided_differences[i][j] << " ";
        //     }
        //     std::cout << "\n";
        // }
    }

  private:
    void computeDividedDifferences() {
        int n = x_values.size();
        for (int i = 0; i < n; ++i) {
            divided_differences[i][0] = f_values[i];
        }

        for (int j = 1; j < n; ++j) {
            for (int i = j; i < n; ++i) {
                divided_differences[i][j] = (divided_differences[i][j - 1] - divided_differences[i - 1][j - 1]) / (x_values[i] - x_values[i - j]);
            }
        }
    }

    std::vector<double> x_values;
    std::vector<double> f_values;
    std::vector<std::vector<double>> divided_differences;
};

class HermiteInterpolator : public Interpolator {
  public:
    HermiteInterpolator(const std::vector<double>& x_values, const std::vector<double>& y_values, const std::vector<double>& yDerivative_values)
        : x_values(x_values), y_values(y_values), yDerivative_values(yDerivative_values), divided_differences(2 * x_values.size(), std::vector<double>(2 * x_values.size(), 0.0)) {
        if (x_values.empty() || y_values.empty() || yDerivative_values.empty() || x_values.size() != y_values.size() || x_values.size() != yDerivative_values.size()) {
            throw std::invalid_argument("x_values, y_values, and yDerivative_values must be non-empty and of the same size");
        }
        computeCoefficients(x_values, y_values, yDerivative_values);
    }

    double interpolate(double x_value) const override {
        int n = x_values.size() * 2;
        double result = divided_differences[0][n - 1];
        for (int i = n - 2; i >= 0; --i) {
            result = divided_differences[i][i] + (x_value - z[i]) * result;
        }
        return result;
    }

    void printInternalData() const {
        std::cout << "x-values: ";
        for (const auto& val : x_values) {
            std::cout << val << " ";
        }
        std::cout << "\nDivided differences: ";

        for (int i = 0; i < divided_differences.size(); ++i) {
            std::cout << divided_differences[i][i] << " ";
        }
        std::cout << "\n";

        // the following code is for debugging purposes
        // std::cout << "\nDivided differences table: \n";
        // for (int i = 0; i < divided_differences.size(); ++i) {
        //     for (int j = 0; j < i + 1; ++j) {
        //         std::cout << divided_differences[i][j] << " ";
        //     }
        //     std::cout << "\n";
        // }
    }

  private:
    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<double> yDerivative_values;
    std::vector<double> z; // Doubled x-values for Hermite interpolation
    std::vector<std::vector<double>> divided_differences;

    void computeCoefficients(const std::vector<double>& x_values, const std::vector<double>& y_values, const std::vector<double>& yDerivative_values) {
        int n = x_values.size();
        z.resize(2 * n);

        for (int i = 0; i < n; ++i) {
            z[2 * i] = x_values[i];
            z[2 * i + 1] = x_values[i];
            divided_differences[2 * i][0] = y_values[i];
            divided_differences[2 * i + 1][0] = y_values[i];
            divided_differences[2 * i + 1][1] = yDerivative_values[i];
            if (i > 0) {
                divided_differences[2 * i][1] = (y_values[i] - y_values[i - 1]) / (x_values[i] - x_values[i - 1]);
            }
        }

        for (int j = 2; j < 2 * n; ++j) {
            for (int i = j; i < 2 * n; ++i) {
                divided_differences[i][j] = (divided_differences[i][j - 1] - divided_differences[i - 1][j - 1]) / (z[i] - z[i - j]);
            }
        }
    }
};

class BezierInterpolator : public PointInterpolator {
  public:
    BezierInterpolator(const std::vector<std::vector<point>>& control_points_sets) : control_points_sets(control_points_sets) {
        if (control_points_sets.empty()) {
            throw std::invalid_argument("Control points sets cannot be empty");
        }
        for (const auto& set : control_points_sets) {
            if (set.size() != 4) {
                throw std::invalid_argument("Each set of control points must contain exactly 4 points for cubic Bezier interpolation");
            }
        }
    }

    point Interpolate(double t) const override {
        if (t < 0.0 || t > 1.0) {
            throw std::invalid_argument("Parameter t must be in the range [0, 1]");
        }

        // Find the appropriate segment for the given t
        double segment_length = 1.0 / control_points_sets.size();
        int segment_index = static_cast<int>(t / segment_length);
        if (segment_index >= control_points_sets.size()) {
            segment_index = control_points_sets.size() - 1;
        }

        // Map t to the range [0, 1]
        double local_t = (t - segment_index * segment_length) / segment_length;

        return calculateBezierPoint(control_points_sets[segment_index], local_t);
    }

  private:
    std::vector<std::vector<point>> control_points_sets;

    point calculateBezierPoint(const std::vector<point>& control_points, double t) const {
        point p0 = control_points[0];
        point p1 = control_points[1];
        point p2 = control_points[2];
        point p3 = control_points[3];

        double u = 1 - t;
        point point_result = p0 * (u * u * u) + p1 * (3 * u * u * t) + p2 * (3 * u * t * t) + p3 * (t * t * t);

        return point_result;
    }
};

#endif