from math import log, log2
from scipy import integrate
import numpy as np


def f(x):
    if x == 0:
        return np.float64(np.pi)
    if x == 1:
        return np.float64(5 * np.pi)
    return np.sin(np.pi * x**5) / (x**5 * (1 - x))


def g(t):
    def _g(x):
        return np.e ** (-np.sqrt(x) + np.sin(x / 10))

    if t + 1e-8 >= 1:
        return 0
    j = 1 / (1 - t) ** 2
    return _g(t / (1 - t)) * j


def simpson(xs, ys):
    N = len(ys) - 1
    h = (xs[-1] - xs[0]) / N
    I1 = h / 3 * sum(ys[i - 1] + 4 * ys[i] + ys[i + 1] for i in range(1, N, 2))
    return I1


def exact_I1():
    N = 1_000_000
    h = 1 / N
    xs = np.array([i * h for i in range(N + 1)])
    ys = list(map(f, xs))
    return integrate.simpson(ys, x=xs)


def exact_I2():
    # taken from wolfram
    return 2.981003452558338


def get_data(a, b, f, N):
    h = (b - a) / N
    xs = np.array([a + i * h for i in range(N + 1)])
    ys = list(map(f, xs))
    return (xs, ys)


def demo_simpson(func, _exact, a, b, e):
    exact = _exact()
    for N in range(20, 50000, 2):
        xs, ys = get_data(a, b, func, N)
        actual = simpson(xs, ys)
        if abs(exact - actual) < e:
            actual_h_over_2 = simpson(*get_data(a, b, func, 2 * N))
            actual_h_over_4 = simpson(*get_data(a, b, func, 4 * N))

            error_h = abs(actual - actual_h_over_2)
            error_h_over_2 = abs(actual_h_over_2 - actual_h_over_4)
            # maybe because of this ratio
            order = log2(error_h / error_h_over_2)
            print(actual)
            print(f"Order of accuracy: {order}")
            return
    print("Something went wrong")


demo_simpson(f, exact_I1, 0, 1, 1e-8)
demo_simpson(g, exact_I2, 0, 1, 1e-8)
