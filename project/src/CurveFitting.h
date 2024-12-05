// CurveFitting.h
// This file defines classes for representing and fitting curves in 2D space.
// Key components:
// - Point: Represents a 2D point with x, y coordinates and provides a method to calculate its length.
// - Curve: Abstract base class for parametric curves, supporting:
//      - operator()(double t): Evaluates the curve at parameter t.
//      - tangent(double t): Computes the tangent vector at t.
//      - plot(): Outputs sampled points of the curve for visualization.
// - CurveFitting: Fits a parametric curve using B-splines for smooth interpolation:
//      - Constructors take a base curve and parameters for fitting.
//      - operator()(double t): Evaluates the fitted curve at parameter t.
//      - Automatically fits Natural or Periodic cubic B-splines for x and y coordinates.

#ifndef CURVEFITTING_H
#define CURVEFITTING_H

#include "BSpline.h"
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

// Abstract curve base class
class Point {
  public:
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    double length() const {
        return std::sqrt(x * x + y * y);
    }
};

const double Pi = 3.14159265358979323846;

// Abstract curve base class
class Curve {
  public:
    double tStart, tEnd;

    Curve(double start = 0, double end = 0) : tStart(start), tEnd(end) {}

    // Pure virtual function for evaluating the curve at parameter t
    virtual Point operator()(double t) const = 0;

    Point tangent(double t) const {
        double h = 1e-6;
        Point p1 = (*this)(t - h);
        Point p2 = (*this)(t + h);
        return {(p2.x - p1.x) / (2 * h), (p2.y - p1.y) / (2 * h)};
    }

    // Function to plot the curve
    void plot(std::ostream& os = std::cout, size_t samples = 300) const {
        if (tEnd <= tStart || samples == 0) {
            throw std::runtime_error("Invalid range or sample count for plotting.");
        }
        double step = (tEnd - tStart) / samples;

        for (size_t i = 0; i < samples; ++i) {
            double t = tStart + i * step;
            Point p = (*this)(t);
            os << p.x << " " << p.y << std::endl;
        }
    }
};

// Curve fitting using B-splines for a parametric 2D curve
class CurveFitting : public Curve {
  public:
    CurveFitting() = default;
    // Constructor that takes a base curve and fitting parameters
    CurveFitting(const Curve& curve, int numControl = 100) : Curve(curve.tStart, curve.tEnd) {
        double speed = (curve.tEnd - curve.tStart) / numControl;
        // Initialize points and cumulative lengths
        double cumulativeLength = 0;
        points_.clear();
        cumulativeLengths_.clear();

        for (double t = curve.tStart; t <= curve.tEnd; t += speed) {
            Point point = curve(t);

            if (!points_.empty()) {
                double dx = point.x - points_.back().x;
                double dy = point.y - points_.back().y;
                cumulativeLength += std::sqrt(dx * dx + dy * dy);
            }

            points_.push_back(point);
            cumulativeLengths_.push_back(cumulativeLength);
        }
        // Ensure proper spline fitting
        fitSplines();
        tStart = 0;
        tEnd = cumulativeLengths_.back();
    }

    CurveFitting(const Curve& curve, double speed = 0.1, bool isEquallySpaced = true) : Curve(curve.tStart, curve.tEnd) {
        if (speed <= 0) {
            throw std::runtime_error("Invalid speed parameter.");
        }
        if (!isEquallySpaced) {
            return;
        }
        // Initialize points and cumulative lengths
        double cumulativeLength = 0;
        points_.clear();
        cumulativeLengths_.clear();
        double t = curve.tStart;
        while (t < curve.tEnd) {
            Point point = curve(t);

            if (!points_.empty()) {
                double dx = point.x - points_.back().x;
                double dy = point.y - points_.back().y;
                cumulativeLength += std::sqrt(dx * dx + dy * dy);
            }
            Point _tangent = curve.tangent(t);
            points_.push_back(point);
            cumulativeLengths_.push_back(cumulativeLength);
            t += speed / _tangent.length();
        }

        t = curve.tEnd;
        Point point = curve(t);
        double dx = point.x - points_.back().x;
        double dy = point.y - points_.back().y;
        cumulativeLength += std::sqrt(dx * dx + dy * dy);
        points_.push_back(point);
        cumulativeLengths_.push_back(cumulativeLength);
        fitSplines();
        tStart = 0;
        tEnd = cumulativeLengths_.back();
    }

    Point operator()(double t) const override {
        if (t < tStart || t > tEnd) {
            throw std::runtime_error("Parameter t out of bounds.");
        }
        double x = splineX_->evaluate(t);
        double y = splineY_->evaluate(t);
        return Point(x, y);
    }

    int control_points_size() const {
        return points_.size();
    }

    // for debugging and checking
    void print_interval_lengths() const {
        for (size_t i = 1; i < cumulativeLengths_.size(); ++i) {
            std::cout << cumulativeLengths_[i] - cumulativeLengths_[i - 1] << std::endl;
        }
    }

    double length() const {
        return tEnd;
    }

  private:
    std::vector<Point> points_;             // Points along the curve
    std::vector<double> cumulativeLengths_; // Cumulative lengths for parameterization
    std::unique_ptr<CubicBSpline> splineX_; // Polymorphic spline for x-coordinate
    std::unique_ptr<CubicBSpline> splineY_; // Polymorphic spline for y-coordinate

    void fitSplines() {
        if (points_.size() < 2) {
            throw std::runtime_error("Insufficient points for spline fitting.");
        }

        // Extract x and y values
        std::vector<double> xValues, yValues;
        for (const Point& point : points_) {
            xValues.push_back(point.x);
            yValues.push_back(point.y);
        }

        // Determine whether to use PeriodicCubicBSpline or NaturalCubicSpline
        bool xPeriodic = std::abs(xValues.front() - xValues.back()) < 1e-6;
        bool yPeriodic = std::abs(yValues.front() - yValues.back()) < 1e-6;

        if (xPeriodic) {
            splineX_ = std::make_unique<PeriodicCubicBSpline>(cumulativeLengths_, xValues);
        } else {
            splineX_ = std::make_unique<NaturalCubicBSpline>(cumulativeLengths_, xValues);
        }

        if (yPeriodic) {
            splineY_ = std::make_unique<PeriodicCubicBSpline>(cumulativeLengths_, yValues);
        } else {
            splineY_ = std::make_unique<NaturalCubicBSpline>(cumulativeLengths_, yValues);
        }
    }
};

class Point3D {
  public:
    double x, y, z;
    Point3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }
};

class SphereCurve {
  public:
    double tStart, tEnd;

    SphereCurve(double start = 0, double end = 0) : tStart(start), tEnd(end) {}

    virtual Point3D operator()(double t) const = 0;

    Point3D tangent(double t) const {
        double h = 1e-6;
        Point3D p1 = (*this)(t - h);
        Point3D p2 = (*this)(t + h);
        return {(p2.x - p1.x) / (2 * h), (p2.y - p1.y) / (2 * h), (p2.z - p1.z) / (2 * h)};
    }

    void plot(std::ostream& os = std::cout, size_t samples = 300) const {
        if (tEnd <= tStart || samples == 0) {
            throw std::runtime_error("Invalid range or sample count for plotting.");
        }
        double step = (tEnd - tStart) / samples;

        for (size_t i = 0; i < samples; ++i) {
            double t = tStart + i * step;
            Point3D p = (*this)(t);
            os << p.x << " " << p.y << " " << p.z << std::endl;
        }
    }
};

class SphereCurvefit : public SphereCurve {
  public:
    SphereCurvefit(const SphereCurve& spCurve, Point3D northPole = {0, 0, 1}, int numControl = 100) : SphereCurve(spCurve.tStart, spCurve.tEnd) {
        build_rotate_matrix(northPole);

        double speed = (spCurve.tEnd - spCurve.tStart) / numControl;
        for (double t = spCurve.tStart; t <= spCurve.tEnd; t += speed) {
            Point3D point = spCurve(t);
            points_.push_back(point);
        }

        fitSplines();
    }

    Point3D operator()(double t) const override {
        if (t < tStart || t > tEnd) {
            throw std::runtime_error("Parameter t out of bounds.");
        }
        Point point = {splineX_->evaluate(t), splineY_->evaluate(t)};
        return inverse_rotate(mapping_to_sphere(point));
    }

  private:
    std::vector<Point3D> points_;
    std::vector<double> cumulativeLengths_;
    std::unique_ptr<CubicBSpline> splineX_;
    std::unique_ptr<CubicBSpline> splineY_;
    std::vector<std::vector<double>> rotate_matrix_ = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector<std::vector<double>> inverse_rotate_matrix_ = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    void build_rotate_matrix(Point3D northPole) {
        double x = northPole.x;
        double y = northPole.y;
        double z = northPole.z;
        double length = northPole.length();

        if (std::abs(length - 1) > 0.001) {
            throw std::runtime_error("Invalid north pole.");
        }

        if (std::abs(x) < 1e-6 && std::abs(y) < 1e-6) {
            if (z > 0) {
                return;
            } else {
                rotate_matrix_ = {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};
                inverse_rotate_matrix_ = {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};
                return;
            }
        }

        double cos_phi = z;
        double sin_phi = std::sqrt(1 - cos_phi * cos_phi);
        double cos_theta = x / sin_phi;
        double sin_theta = y / sin_phi;

        rotate_matrix_ = {
            {cos_theta * cos_phi, sin_theta * cos_phi, -sin_phi},
            {-sin_theta, cos_theta, 0},
            {cos_theta * sin_phi, sin_theta * sin_phi, cos_phi},
        };

        inverse_rotate_matrix_ = {
            {cos_theta * cos_phi, -sin_theta, cos_theta * sin_phi},
            {sin_theta * cos_phi, cos_theta, sin_theta * sin_phi},
            {-sin_phi, 0, cos_phi}};
    }

    Point3D rotate(Point3D point) {
        double x = point.x;
        double y = point.y;
        double z = point.z;
        return {rotate_matrix_[0][0] * x + rotate_matrix_[0][1] * y + rotate_matrix_[0][2] * z,
                rotate_matrix_[1][0] * x + rotate_matrix_[1][1] * y + rotate_matrix_[1][2] * z,
                rotate_matrix_[2][0] * x + rotate_matrix_[2][1] * y + rotate_matrix_[2][2] * z};
    }

    Point3D inverse_rotate(Point3D point) const {
        double x = point.x;
        double y = point.y;
        double z = point.z;
        return {inverse_rotate_matrix_[0][0] * x + inverse_rotate_matrix_[0][1] * y + inverse_rotate_matrix_[0][2] * z,
                inverse_rotate_matrix_[1][0] * x + inverse_rotate_matrix_[1][1] * y + inverse_rotate_matrix_[1][2] * z,
                inverse_rotate_matrix_[2][0] * x + inverse_rotate_matrix_[2][1] * y + inverse_rotate_matrix_[2][2] * z};
    }

    Point mapping_to_plain(Point3D point) const {
        double x = point.x;
        double y = point.y;
        double z = point.z;
        return {x / (1 - z), y / (1 - z)};
    }

    Point3D mapping_to_sphere(Point point) const {
        double x = point.x;
        double y = point.y;
        double denominator = 1 + x * x + y * y;
        return {2 * x / denominator, 2 * y / denominator, (-1 + x * x + y * y) / denominator};
    }

    void fitSplines() {
        if (points_.size() < 2) {
            throw std::runtime_error("Insufficient points for spline fitting.");
        }

        std::vector<Point> mapped_control_points;
        for (const Point3D& point : points_) {
            Point3D rotatedPoint = rotate(point);
            Point mappedPoint = mapping_to_plain(rotatedPoint);
            mapped_control_points.push_back(mappedPoint);
        }

        // Extract x and y values
        std::vector<double> xValues, yValues;
        double cumulativeLength = 0;
        for (const Point& point : mapped_control_points) {
            if (!points_.empty()) {
                double dx = point.x - points_.back().x;
                double dy = point.y - points_.back().y;
                cumulativeLength += std::sqrt(dx * dx + dy * dy);
            }
            xValues.push_back(point.x);
            yValues.push_back(point.y);
            cumulativeLengths_.push_back(cumulativeLength);
        }

        // Determine whether to use PeriodicCubicBSpline or NaturalCubicSpline
        bool xPeriodic = std::abs(xValues.front() - xValues.back()) < 1e-6;
        bool yPeriodic = std::abs(yValues.front() - yValues.back()) < 1e-6;

        if (xPeriodic) {
            splineX_ = std::make_unique<PeriodicCubicBSpline>(cumulativeLengths_, xValues);
        } else {
            splineX_ = std::make_unique<NaturalCubicBSpline>(cumulativeLengths_, xValues);
        }

        if (yPeriodic) {
            splineY_ = std::make_unique<PeriodicCubicBSpline>(cumulativeLengths_, yValues);
        } else {
            splineY_ = std::make_unique<NaturalCubicBSpline>(cumulativeLengths_, yValues);
        }

        tStart = 0;
        tEnd = cumulativeLengths_.back();
    }
};

#endif // CURVEFITTING_H
