#ifndef CURVEFITTING_H
#define CURVEFITTING_H

#include "BSpline.h"
#include <cmath>
#include <stdexcept>
#include <vector>

// Representation of a 2D point
class Point {
  public:
    double x, y;
    Point(double x, double y) : x(x), y(y) {}
};

// Abstract curve base class
class Curve {
  public:
    virtual Point operator()(double t) const = 0;
};

const double Pi = 3.14159265358979323846;

// Curve fitting using B-splines for a parametric 2D curve
class CurveFitting : public Curve {
  public:
    // Constructor that takes a vector of points
    CurveFitting(const Curve& curve, double speed = 0.1, double start = 0, double end = 2 * Pi) {
        // Compute chordal lengths
        double cumulativeLength = 0;
        for (double t = start; t < end; t += speed) {
            Point point = curve(t);
            double dx = point.x - points_.back().x;
            double dy = point.y - points_.back().y;
            cumulativeLength += std::sqrt(dx * dx + dy * dy);
            points_.push_back(point);
            cumulativeLengths_.push_back(cumulativeLength);
        }

        // Fit splines for x and y coordinates
        fitSplines();
    }

  private:
    std::vector<Point> points_;             // points
    std::vector<double> cumulativeLengths_; // Cumulative chordal lengths
    PeriodicCubicBSpline splineX_;          // Spline for x-coordinate
    PeriodicCubicBSpline splineY_;          // Spline for y-coordinate

    // Fit splines for the x and y coordinates
    void fitSplines() {
        // Extract x and y values from the points
        std::vector<double> xValues;
        std::vector<double> yValues;
        for (const Point& point : points_) {
            xValues.push_back(point.x);
            yValues.push_back(point.y);
        }

        xValues.push_back(points_.front().x);
        yValues.push_back(points_.front().y);

        // Create cubic B-splines for x and y coordinates
        splineX_ = PeriodicCubicBSpline(cumulativeLengths_, xValues);
        splineY_ = PeriodicCubicBSpline(cumulativeLengths_, yValues);
    }
};

#endif // CURVEFITTING_H
