#include "Interpolator.hpp"
#include <cmath>
#include <fstream>
#include <vector>

struct pointWithTangent {
    point position;
    point tangent;

    pointWithTangent(const point& pos, const point& tan) : position(pos), tangent(tan) {}
};

class Heartfunction {
  public:
    point operator()(double t) const {
        double x = std::sqrt(3) * std::sin(M_PI * t);
        double y = (2.0 / 3.0) * (std::sqrt(3) * std::cos(M_PI * t) + std::sqrt(x));
        return point(x, y);
    }

    point tangent(double t) const {
        double x = std::sqrt(3) * std::sin(M_PI * t);
        double dx_dt = M_PI * std::sqrt(3) * std::cos(M_PI * t);
        double dy_dt = (2.0 / 3.0) * (-M_PI * std::sqrt(3) * std::sin(M_PI * t) + (1.0 / (2.0 * std::sqrt(x) + 1e-8)) * dx_dt);
        point tangent(dx_dt, dy_dt);
        tangent.normalize();
        return tangent;
    }
};

std::vector<pointWithTangent> computeMarkerPointsWithTangents(int m) {
    std::vector<pointWithTangent> marker_points_with_tangents;
    double t_step = 1.0 / m;

    for (int i = 0; i <= m; ++i) {
        double t = i * t_step;

        Heartfunction heartfunction;

        point position = heartfunction(t);
        point tangent = heartfunction.tangent(t);

        marker_points_with_tangents.emplace_back(position, tangent);
    }

    return marker_points_with_tangents;
}

std::vector<std::vector<point>> computeControlPoints(const std::vector<pointWithTangent>& marker_points_with_tangents) {
    std::vector<std::vector<point>> control_points_sets;

    for (int j = 0; j < marker_points_with_tangents.size() - 1; ++j) {
        point p_j = marker_points_with_tangents[j].position;
        point p_j1 = marker_points_with_tangents[j + 1].position;

        // set the tangent vectors so that their lengths are equal to the distance between the points
        point tangent_j = marker_points_with_tangents[j].tangent;
        point tangent_j1 = marker_points_with_tangents[j + 1].tangent;
        tangent_j = tangent_j * distance(p_j, p_j1);
        tangent_j1 = tangent_j1 * distance(p_j, p_j1);

        point q_0 = p_j;
        point q_1 = (tangent_j * (1.0 / 3.0)) + p_j;
        point q_2 = p_j1 - (tangent_j1 * (1.0 / 3.0));
        point q_3 = p_j1;

        control_points_sets.push_back({q_0, q_1, q_2, q_3});
    }

    return control_points_sets;
}

int main() {
    std::ofstream outFile("plotCode/pointsInF.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error opening file for writing.\n";
        return 1;
    }

    Heartfunction heartfunction;

    std::vector<int> m_values = {5, 20, 80};
    for (int m : m_values) {
        outFile << "For m = " << m << ":\n";

        std::vector<pointWithTangent> marker_points_with_tangents = computeMarkerPointsWithTangents(m);
        std::vector<std::vector<point>> control_points_sets = computeControlPoints(marker_points_with_tangents);
        BezierInterpolator bezier_interpolator(control_points_sets);

        double max_distance = 0.0;
        double max_distance_for_middle_t = 0.0;
        for (double t = 0; t <= 1.0; t += 0.001) {
            point interpolated_point = bezier_interpolator.Interpolate(t);
            outFile << interpolated_point.getX() << " " << interpolated_point.getY() << "\n"; // Write x, y coordinates

            // compute the distance between the interpolated point and the actual point
            double dis = distance(interpolated_point, heartfunction(t));
            if (dis > max_distance) {
                max_distance = dis;
            }
            if (t > 0.2 && t < 0.8) {
                if (dis > max_distance_for_middle_t) {
                    max_distance_for_middle_t = dis;
                }
            }
        }

        outFile << "\n";

        std::cout << "For m = " << 2 * m << ", The maximum distance between the interpolated points and the actual points is " << max_distance << std::endl;
        std::cout << "For m = " << 2 * m << ", The maximum distance between the interpolated points and the actual points for t in (0.2, 0.8) is " << max_distance_for_middle_t << std::endl;
    }

    outFile.close();
    std::cout << "Data written to file.\n";
    std::cout << "Run the following command in the terminal to plot the data:\n";
    std::cout << "python3 plotCode/plotF.py (you can replace python3 with the Python interpreter on your computer)\n";

    return 0;
}
