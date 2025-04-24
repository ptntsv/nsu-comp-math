from scipy.interpolate import CubicSpline
import numpy as np


def f(x):
    return 1 / (1 + 25 * x**2)


def split_segment(l, r, n):
    step = (r - l) / n
    return [l + i * step for i in range(n)]


def setup1():
    rel = {}
    for n in range(7, 20, 2):
        x = [2 * i / n - 1 for i in range(n + 1)]
        y = list(map(f, x))
        cs = CubicSpline(x, y, bc_type=((2, 0), (2, 0)))
        k = 10
        xs = list()
        for i in range(1, len(x)):
            xs += split_segment(x[i - 1], x[i], k)
        print(cs(xs))
        return
        # xs.append(x[-1])
        # actual = cs(xs)
        # expected = list(map(f, xs))
        # error = np.average([abs(actual[i] - expected[i]) for i in range(len(actual))])
        # rel[n] = error
    print(rel)


def setup2():
    x = [0, 1, 2, 3]
    y = [0, 1, 0, -1]
    cs = CubicSpline(x, y, bc_type=((2, 1), (2, 2)))
    xs = np.arange(0, 3.1, 0.1)
    print(cs(xs))


setup1()
