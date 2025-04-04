import matplotlib.pyplot as plt
import numpy as np
import math


class newton_polynomial:
    def __init__(self, xs, ys):
        self.xs = xs
        self.cfs = [ys[0]] + self._poly_newton_coefficient(list(zip(xs, ys)))

    def eval(self, x):
        ans = 0
        for i in range(len(self.cfs) - 1, -1, -1):
            ans = self.cfs[i] + ans * (x - self.xs[i])
        return ans

    def _poly_newton_coefficient(self, points):
        n = len(points)

        def div_diff(x1, y1, x0, y0):
            return (y1 - y0) / (x1 - x0)

        dp = [[0] * n for _ in range(n)]
        for k in range(1, n):
            for i in range(n - k):
                if k == 1:
                    dp[i][i + k] = div_diff(*points[i + k], *points[i])
                else:
                    j = i + k
                    dp[i][j] = (dp[i + 1][j] - dp[i][j - 1]) / (
                        points[j][0] - points[i][0]
                    )
        return dp[0][1:]


def f(x):
    return 1 / (1 + 25 * x**2)


def interpolate(f, data, xs):
    ydata = list(map(f, data))
    p = newton_polynomial(data, ydata)
    return list(map(p.eval, xs))


def demo(a=-1, b=1, f=f):
    x0 = -1.0
    x1 = 1.0
    step = 0.001
    _, axs = plt.subplots(3)

    axs[0].set_title("Default")
    axs[1].set_title("Chebyshev polynomial")

    axs[0].label_outer()
    axs[1].label_outer()
    for n in range(3, 10, 2):
        xs = np.arange(x0, x1, step)
        ys = list(map(f, xs))
        chebyshev_data = [
            0.5 * (b + a)
            + 0.5 * (b - a) * math.cos((2 * k + 1) * math.pi / (2 * n + 2))
            for k in range(n + 1)
        ]

        xdata = [2 * i / n - 1 for i in range(n + 1)]
        print(xdata)

        pvalues_default = interpolate(f, xdata, xs)
        pvalues_chebyshev = interpolate(f, chebyshev_data, xs)

        axs[0].plot(xs, [abs(pvalues_default[i] - ys[i]) for i in range(len(ys))])
        axs[1].plot(xs, [abs(pvalues_chebyshev[i] - ys[i]) for i in range(len(ys))])
    plt.show()


demo()
