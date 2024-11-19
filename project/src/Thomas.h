#ifndef THOMAS_H
#define THOMAS_H

#include <cmath>
#include <vector>

// a: sub-diagonal, b: main diagonal, c: super-diagonal, d: right-hand side
void thomasAlgorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x) {
    // cope with the case where b[0] is zero, which is needed in cyclicthomasAlgorithm
    if (std::abs(b[0]) < 1e-8) {
        std::vector<double> d_star(d.begin() + 1, d.end());
        double x2 = d[0] / c[0];
        d_star[0] -= x2 * b[1];
        d_star[1] -= x2 * a[1];
        std::vector<double> a_star(a.begin() + 1, a.end());
        a_star[0] = 0;
        std::vector<double> b_star(b.begin() + 1, b.end());
        b_star[0] = a[0];
        std::vector<double> c_star(c.begin() + 1, c.end());
        thomasAlgorithm(a_star, b_star, c_star, d_star, x);
        x.insert(x.begin() + 1, x2);
        x.pop_back();
        return;
    }
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

void cyclicthomasAlgorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x) {
    int n = b.size();
    std::vector<double> u(n, 0.0);
    std::vector<double> v(n, 0.0);
    u[0] = b[0];
    v[0] = 1;
    u[n - 1] = c[n - 1];
    v[n - 1] = a[0] / b[0];
    std::vector<double> b_star = b;
    b_star[n - 1] -= a[0] * c[n - 1] / b_star[0];
    b_star[0] = 0;
    std::vector<double> y(n, 0.0);
    std::vector<double> q(n, 0.0);
    std::vector<double> a_star(a.begin() + 1, a.end());
    std::vector<double> c_star(c.begin(), c.end() - 1);
    thomasAlgorithm(a_star, b_star, c_star, d, y);
    thomasAlgorithm(a_star, b_star, c_star, u, q);
    double z = (v[0] * y[0] + v[n - 1] * y[n - 1]) / (1 + q[0] * v[0] + q[n - 1] * v[n - 1]);
    for (int i = 0; i < n; ++i) {
        x[i] = y[i] - z * q[i];
    }
}

#endif // THOMAS_H
