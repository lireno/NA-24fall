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

        for (size_t i = 0; i <= samples; ++i) {
            double t = tStart + i * step;
            Point p = (*this)(t);
            os << p.x << " " << p.y << std::endl;
        }
    }
};

// Curve fitting using B-splines for a parametric 2D curve
class CurveFitting : public Curve {
  public:
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
    }

    Point operator()(double t) const override {
        if (t < tStart || t > tEnd) {
            throw std::runtime_error("Parameter t out of bounds.");
        }
        double x = splineX_->evaluate(t);
        double y = splineY_->evaluate(t);
        return Point(x, y);
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

#endif // CURVEFITTING_H
