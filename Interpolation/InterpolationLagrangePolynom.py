import numpy as np
import sympy as sp
import mpmath as mp
import scipy
import matplotlib.pyplot as plt

n = 11
a, b = 1, 2
xSymbol = sp.S('x')
f = lambda x: 0.8 * np.log(x) + 0.6 * np.cos(x)
x = np.linspace(a, b, n)
y = f(x)

p = 0
Q = [0] * n
for k in range(n):
    dw = np.prod(sp.S('x') - x) / (sp.S('x') - x[k])
    Q[k] = dw / dw.subs('x', x[k])
    p += y[k] * Q[k]

xRange = np.linspace(a - (b - a) / 10, b + (b - a) / 10, (b - a) * 100)
print("P(x) =", sp.expand(p))
print("Accuracy: %e" % np.max([np.abs(f(value) - p.subs('x', value)) for value in xRange]))
plt.plot(x, y, 'ro')
plt.plot(xRange, [p.subs(xSymbol, value) for value in xRange])
plt.show()