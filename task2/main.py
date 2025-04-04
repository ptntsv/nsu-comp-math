import matplotlib.pyplot as plt
import numpy as np


def f(z: complex):
    return z**3 - 1.0


def df(z: complex):
    return 3.0 * z**2


def root_color(root, e=0.0001):
    if abs(root - 1) < e:
        return "b"
    if abs(root - complex(-0.5000, -0.86603)) < e:
        return "g"
    return "r"


def newton(guess, ex, max_iters=1000000) -> complex:
    pguess = guess
    for _ in range(max_iters):
        guess = pguess - f(pguess) / df(pguess)
        if abs(pguess - guess) < ex:
            break
        pguess = guess
    return complex(round(guess.real, 4), round(guess.imag, 4))


def plot_newton_fractal(step=0.01, ex=0.00001):
    side = 1
    re = np.arange(-side, side, step)
    im = np.arange(-side, side, step)
    for x in re:
        for y in im:
            root = newton(complex(x, y), ex)
            plt.plot(x, y, root_color(root) + "o")

    plt.show()


plot_newton_fractal()
