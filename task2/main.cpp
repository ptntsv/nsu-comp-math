#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

using namespace std;
using comp = std::complex<double>;

double f(const double x) { return tan(x) - x; }
double df(const double x) { return 1 / (cos(x) * cos(x)) - 1; }

comp g(const comp z) { return z * z * z - 1.0; }
comp dg(const comp z) { return z * z * 3.0; }

double bisection(double left, double right, double accuracy) {
    double a, b, c = 0;
    a = left, b = right;
    double fa = f(a);
    while (std::abs(b - a) > 2 * accuracy) {
        c = (a + b) / 2.0;
        double fc = f(c);
        if (fc == 0)
            break;
        if (fa * fc < 0) {
            b = c;
        } else {
            a = c;
            fa = fc;
        }
    }
    return c;
}

double fixed_point(double l, double r, double e, int period,
                   int maxIters = 10000) {
    double x = (l + r) / 2;
    double atanx;
    for (size_t i = 0; i < maxIters; i++, x = atanx) {
        atanx = atan(x) + M_PI * period;
        if (std::abs(atanx - x) < e)
            return atanx;
    }
    return x;
}

template <typename T>
T newton(T guess, double ex, T (*f)(const T), T (*df)(const T),
         size_t maxIters = 1000000) {
    T pguess = guess;
    for (size_t i = 0; i < maxIters; pguess = guess, i++) {
        guess = pguess - f(pguess) / df(pguess);
        if (std::abs(pguess - guess) < ex)
            break;
    }
    return guess;
}

void secant_step(double& a, double& b, double& delta) {
    double tmp = b;
    b = b - f(b) * (b - a) / (f(b) - f(a));
    a = tmp;
    delta = std::abs(a - b);
}

double secant(double a, double b, double e) {
    double delta;
    do {
        secant_step(a, b, delta);
    } while (delta > e);

    double prev_delta = delta;
    do {
        secant_step(a, b, delta);
    } while (delta < prev_delta && delta > 0);

    return b;
}

void demo() {
    int periods = 10;
    double epsilon = 0.001;
    // cout << "===============" << endl;
    // cout << "Bisection method:" << endl;
    // for (int i = 0; i < periods; i++) {
    //     cout << bisection(M_PI * i - M_PI_2 + epsilon,
    //                       M_PI * i + M_PI_2 - epsilon, 0.0001)
    //          << endl;
    // }
    cout << "===============" << endl;
    cout << "Fixed point iteration method" << endl;
    for (int i = 0; i < periods; i++) {
        cout << fixed_point(M_PI * i - M_PI_2 + epsilon,
                            M_PI * i + M_PI_2 - epsilon, 0.0001, i)
             << endl;
    }
    // cout << "===============" << endl;
    // cout << "Newton's method" << endl;
    // cout << newton<double>(7.8, 0.00001, &f, &df) << endl;
    //
    // cout << "===============" << endl;
    // cout << "Newton's method(complex)" << endl;
    // cout << newton<comp>(comp{2, 1}, 0.00001, &g, &dg) << endl;
    //
    // cout << "===============" << endl;
    // cout << "Secant method" << endl;
    // cout << secant(M_PI_2 + 0.5, 4.5, 0.0001);
}

int main() {
    demo();
    return 0;
}
