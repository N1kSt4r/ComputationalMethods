import numpy as np
import scipy as sp
from scipy.misc import derivative
import matplotlib.pyplot as plt

n = 5
a, b = 1, 2
f = lambda x: 0.8 * np.log(x) + 0.6 * np.cos(x)
x = np.linspace(a, b, n)
y = f(x)

w = lambda arg: np.prod(arg - x)
dw = lambda arg, k: np.prod(arg - np.append(x[:k], x[k+1:]))

def P(arg):
    value = 0
    for k in range(len(x)):
        value += y[k] * dw(arg, k) / dw(x[k], k)
    return value

def R(arg):
    max_derivative = max([
        abs(derivative(f, n=n+1, x0=1, dx=0.1, order=n//2 * 2 + 5))
        for i in np.linspace(a, b, n * 100)
    ])
    return abs(w(arg) * max_derivative / np.prod(np.arange(n + 1) + 1))

xRange = np.linspace(a - (b - a) / 10, b + (b - a) / 10, (b - a) * 100)
print("Value at (b + a)/5: %f" % P((b + a) / 5))
print("Accuracy at (b + a)/5: %e" % R((b + a) / 5))
plt.plot(x, y, 'ro')
plt.plot(xRange, [P(value) for value in xRange])
plt.show()