from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt

_, axs = plt.subplots(2)
axs[0].set_title("Actual")
axs[1].set_title("Expected")


def f(x):
    return 1 / (1 + 25 * x**2)


def split_segment(l, r, n):
    step = (r - l) / n
    return [l + i * step for i in range(n)]


def setup1():
    rel = {}
    for n in range(5, 13, 2):
        x = [2 * i / n - 1 for i in range(n + 1)]
        y = list(map(f, x))
        cs = CubicSpline(x, y, bc_type=((2, 0), (2, 0)))
        k = 50
        xs = list()
        for i in range(1, len(x)):
            xs += split_segment(x[i - 1], x[i], k)
        xs.append(x[-1])
        actual = cs(xs)
        expected = list(map(f, xs))
        axs[1].plot(xs, [abs(actual[i] - expected[i]) for i in range(len(expected))])


def setup2():
    x = [0, 1, 2, 3]
    y = [0, 1, 0, -1]
    cs = CubicSpline(x, y, bc_type=((2, 1), (2, 2)))
    xs = np.arange(0, 3.1, 0.1)
    print(cs(xs))


def plot():
    nmap = dict()
    with open("../result.txt", "r") as file:
        for line in file:
            items = line.strip().split(" ")
            n = int(items[0])
            x = float(items[1])
            diff = float(items[2])

            if n not in nmap:
                nmap[n] = ([], [])
            nmap[n][0].append(x)
            nmap[n][1].append(diff)

    for k, v in nmap.items():
        axs[0].plot(v[0], v[1], label=f"{k}")

    # xs, actual, expected_ys = setup1()
    # print(xs[0], actual[0], expected_ys[0])
    # plt.plot(xs, expected_ys, label="Actual")
    plt.legend(loc="upper right")


setup1()
plot()
plt.show()
