#include <assert.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
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

vector<double> get_spline_coefficients(vector<double>& x, vector<double>& f,
                                       double c0, double cn) {
    int n = x.size() - 1;  // number of segments
    vector<double> h(n);
    int m = n - 1;  // SLE matrix dimension
    for (size_t i = 0; i < n; i++) h[i] = x[i + 1] - x[i];
    vector<double> A(m), B(m), C(m), D(m);
    for (size_t i = 0; i < m; i++) {
        if (i == 0) {
            A[i] = 0;
            B[i] = (h[0] + h[1]) / 3;
            C[i] = h[1] / 6;
            D[i] = (f[2] - f[1]) / h[1] - (f[1] - f[0]) / h[0] - h[0] * c0 / 6;
        } else if (i == m - 1) {
            A[i] = h[i] / 6;
            B[i] = (h[i] + h[i + 1]) / 3;
            C[i] = 0;
            D[i] = (f[i + 2] - f[i + 1]) / h[i + 1] - (f[i + 1] - f[i]) / h[i] -
                   h[i + 1] * cn / 6;
        } else {
            A[i] = h[i] / 6.0;
            B[i] = (h[i] + h[i + 1]) / 3.0;
            C[i] = h[i + 1] / 6.0;
            D[i] = (f[i + 2] - f[i + 1]) / h[i + 1] - (f[i + 1] - f[i]) / h[i];
        }
    }

    vector<double> c = thomas(A, B, C, D);
    vector<double> full(n + 1);
    full[0] = c0;
    full[n] = cn;
    for (int i = 1; i < n; ++i) {
        full[i] = c[i - 1];
    }
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
            A[i] = 0;
            B[i] = h[i] / 3;
            C[i] = h[i] / 6;
            D[i] = (y[1] - y[0]) / h[0] - lCond;
        } else if (i == m - 1) {
            A[i] = h[i - 1] / 6;
            B[i] = h[i - 1] / 3;
            C[i] = 0;
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

void interpolate(vector<double>& x, std::function<double(double)> f, double c0,
                 double cn) {
    vector<double> y(x.size());
    std::transform(x.begin(), x.end(), y.begin(), f);
    map<double, double> interpolated{};
    int m = x.size() - 1;
    int k = 10;  // number of points inside segment
    vector<double> h(m);
    // vector<double> c(get_spline_coefficients(x, y, c0, cn));
    vector<double> c(get_spline_coefficients_2(x, y, c0, cn));
    for (size_t i = 0; i < m; i++) h[i] = x[i + 1] - x[i];
    for (size_t i = 0; i < x.size() - 1; i++) {
        vector<double> subseg = split_segment(x[i], x[i + 1], k);
        for (auto& s : subseg) {
            double spline = eval_spline_polynomial(
                s, make_pair(y[i], y[i + 1]), make_pair(x[i], x[i + 1]), h[i],
                make_pair(c[i], c[i + 1]));
            interpolated.insert(make_pair(s, spline));
        }
    }
    interpolated.insert(make_pair(x[x.size() - 1], y[y.size() - 1]));
    for (auto& tuple : interpolated) {
        cout << setprecision(7) << tuple.second << " ";
    }
    cout << endl;
}

int main() {
    for (int n = 7; n < 10; n += 2) {
        vector<double> x(n + 1);
        for (size_t i = 0; i < n + 1; i++) x[i] = 2.0 * i / n - 1;
        interpolate(x, f, 4.0, 4.5);
        break;
    }
}
