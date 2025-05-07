#include <assert.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

using namespace std;

template <typename T>
void println(const vector<T>& v) {
    for (auto& x : v) {
        cout << setprecision(7) << x << " ";
    }
    cout << endl;
}

inline double f(double x) { return 1 / (1 + 25 * x * x); }

vector<double> thomas(vector<double>& A, vector<double>& B, vector<double>& C,
                      vector<double>& D) {
    int m = A.size();
    vector<double> alpha(m), beta(m);
    alpha[0] = -C[0] / B[0];
    beta[0] = D[0] / B[0];

    for (size_t i = 1; i < m; i++) {
        double y = B[i] + A[i] * alpha[i - 1];
        alpha[i] = -C[i] / y;
        beta[i] = (D[i] - A[i] * beta[i - 1]) / y;
    }

    vector<double> coefficients(m);
    coefficients[m - 1] = beta[m - 1];

    for (int i = m - 2; i >= 0; --i) {
        coefficients[i] = alpha[i] * coefficients[i + 1] + beta[i];
    }

    return coefficients;
}

vector<double> get_spline_coefficients(vector<double>& x, vector<double>& y,
                                       double c0, double cn) {
    int n = x.size() - 1;  // number of segments
    vector<double> h(n);
    int m = n - 1;  // SLE matrix dimension
    for (size_t i = 0; i < n; i++) h[i] = x[i + 1] - x[i];
    vector<double> A(m), B(m), C(m), D(m);

    for (size_t i = 0; i < m; i++) {
        if (i == 0) {
            B[i] = 2 * (h[0] + h[1]);
            C[i] = h[1];
            D[i] =
                6 * ((y[2] - y[1]) / h[1] - (y[1] - y[0]) / h[0]) - h[0] * c0;
        } else if (i == m - 1) {
            A[i] = h[i];
            B[i] = 2 * (h[i] + h[i + 1]);
            D[i] = 6 * ((y[i + 2] - y[i + 1]) / h[i + 1] -
                        (y[i + 1] - y[i]) / h[i]) -
                   h[i + 1] * cn;
        } else {
            A[i] = h[i];
            B[i] = 2 * (h[i] + h[i + 1]);
            C[i] = h[i + 1];
            D[i] = 6 * ((y[i + 2] - y[i + 1]) / h[i + 1] -
                        (y[i + 1] - y[i]) / h[i]);
        }
    }
    vector<double> c = thomas(A, B, C, D);
    vector<double> full(n + 1);
    for (size_t i = 1; i < n; i++) {
        full[i] = c[i - 1];
    }
    full[0] = c0;
    full[n] = cn;
    return full;
}

vector<double> get_spline_coefficients_2(vector<double>& x, vector<double>& y,
                                         double lCond, double rCond) {
    int m = x.size();
    vector<double> h(m - 1);
    for (size_t i = 0; i < m - 1; i++) h[i] = x[i + 1] - x[i];
    vector<double> A(m), B(m), C(m), D(m);
    for (size_t i = 0; i < m; i++) {
        if (i == 0) {
            B[i] = h[i] / 3;
            C[i] = h[i] / 6;
            D[i] = (y[1] - y[0]) / h[0] - lCond;
        } else if (i == m - 1) {
            A[i] = h[i - 1] / 6;
            B[i] = h[i - 1] / 3;
            D[i] = rCond - (y[i] - y[i - 1]) / h[i - 1];
        } else {
            A[i] = h[i - 1] / 6;
            B[i] = (h[i - 1] + h[i]) / 3;
            C[i] = h[i] / 6;
            D[i] = (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1];
        }
    }
    return thomas(A, B, C, D);
}

double eval_spline_polynomial(double p, pair<double, double> y,
                              pair<double, double> x, double h,
                              pair<double, double> c) {
    double k = y.first / h * (x.second - p);
    double l = y.second / h * (p - x.first);
    double m =
        c.first / (6 * h) * (pow(x.second - p, 3) - h * h * (x.second - p));
    double n =
        c.second / (6 * h) * (pow(p - x.first, 3) - h * h * (p - x.first));
    return k + l + m + n;
}

vector<double> split_segment(double begin, double end, int n) {
    vector<double> subsegments(n);
    double step = (end - begin) / n;
    for (int i = 0; i < n; i++) subsegments[i] = begin + i * step;
    return subsegments;
}

map<double, double> interpolate(vector<double>& x,
                                std::function<double(double)> f, double c0,
                                double cn) {
    vector<double> y(x.size());
    std::transform(x.begin(), x.end(), y.begin(), f);
    map<double, double> interpolated{};
    int m = x.size() - 1;
    int k = 50;  // number of points inside segment
    vector<double> h(m);
    vector<double> c(get_spline_coefficients(x, y, c0, cn));
    // vector<double> c(get_spline_coefficients_2(x, y, c0, cn));
    for (size_t i = 0; i < m; i++) h[i] = x[i + 1] - x[i];
    for (size_t i = 0; i < x.size() - 1; i++) {
        vector<double> subseg = split_segment(x[i], x[i + 1], k);
        if (i + 1 == x.size() - 1) subseg.push_back(x[x.size() - 1]);
        for (auto& s : subseg) {
            double spline = eval_spline_polynomial(
                s, make_pair(y[i], y[i + 1]), make_pair(x[i], x[i + 1]), h[i],
                make_pair(c[i], c[i + 1]));
            interpolated.insert(make_pair(s, spline));
        }
    }
    return interpolated;
}

double avg_error(vector<double>& xs, std::function<double(double)> f,
                 double lcond, double rcond) {
    int n = xs.size();
    map<double, double> actual = interpolate(xs, f, lcond, rcond);
    vector<double> errors;
    for (auto const& [x, y] : actual) {
        cout << n - 1 << " " << x << " " << std::abs(y - f(x)) << endl;
        errors.push_back(abs(y - f(x)));
    }
    return accumulate(errors.begin(), errors.end(), 0.0) / errors.size();
}

void cubic_spline_demo() {
    double lcond = 0, rcond = 0;
    map<int, double> error_rel{};
    for (int n = 5; n < 12; n += 2) {
        vector<double> x(n + 1);
        for (size_t i = 0; i < n + 1; i++) x[i] = 2.0 * i / n - 1;
        error_rel.insert(make_pair(n, avg_error(x, f, lcond, rcond)));
    }
    // for (auto const& [n, error] : error_rel) {
    //     cout << n << ": " << error << endl;
    // }
}

void first_row(int order, double value, double y, double h, double& a,
               double& b, double& c, double& d) {
    if (order == 0) {
        b = -2;
        c = 1;
        d = h * h * y - value;
    } else if (order == 1) {
        b = 1;
        c = -1;
        d = -h * (value + h * y);
    }
}

void last_row(int order, double value, double y, double h, double& a, double& b,
              double& c, double& d) {
    if (order == 0) {
    }
}

vector<double> solve(const vector<double>& y, pair<int, int> order,
                     pair<double, double> values, double h) {
    int m = y.size() - 2;
    double lcond, rcond;
    vector<double> A(m), B(m), C(m), D(m);
    for (size_t i = 0; i < m; i++) {
        if (i == 0) {
            first_row(order.second, values.first, y[i], h, A[i], B[i], C[i],
                      D[i]);
        } else if (i == m - 1) {
            last_row(order.second, values.second, y[i], h, A[i], B[i], C[i],
                     D[i]);
            A[i] = 1;
            B[i] = -2;
            C[i] = 1;
            D[i] = h * h * y[i] - C[i] * rcond;
        } else {
            A[i] = 1;
            B[i] = -2;
            C[i] = 1;
            D[i] = h * h * y[i];
        }
    }
    vector<double> roots = thomas(A, B, C, D);

    int n = m + 2;
    vector<double> full(n);
    for (size_t i = 1; i < full.size(); i++) {
        full[i] = roots[i - 1];
    }

    full[0] = lcond;
    full[n - 1] = rcond;

    return full;
}

double u(double x) { return cos(x); }
double u_actual(double x) { return -1 + x * x / 2 - pow(x, 4) / 24; }

void de_demo() {
    double offset = 1e-8;
    pair<int, int> order{0, 0};
    pair<double, double> values{1, 2};
    double k = 10;
    double a = -M_PI_2, b = M_PI_2;
    vector<double> xs = split_segment(a, b, k);
    double h = (b - a - offset) / k;
    vector<double> ys(xs.size());
    transform(xs.begin(), xs.end(), ys.begin(), u);
    vector<double> solution = solve(ys, order, values, h);
    vector<double> u_act(xs.size());
    transform(xs.begin(), xs.end(), u_act.begin(), u_actual);
    println(solution);
    println(u_act);
}

int main() { cubic_spline_demo(); }
