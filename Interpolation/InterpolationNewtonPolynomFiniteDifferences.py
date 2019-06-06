import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

n = 11
a, b = 1, 2
xSymbol = sp.S('x')
f = lambda x: 0.8 * np.log(x) + 0.6 * np.cos(x)
x = np.linspace(a, b, n)
y = f(x)

diff = np.zeros((n, n))
diff[0] = np.copy(y)

for k in range(1, n):
    for i in range(n - k):
        diff[k][i] = (diff[k - 1][i + 1] - diff[k - 1][i]) / (x[i + k] - x[i])

p = 0
w = 1
for i in range(n - 1):
    p += diff[i][0] * w
    w *= xSymbol - x[i]

xRange = np.linspace(a - (b - a) / 10, b + (b - a) / 10, (b - a) * 100)
print("P(x) =", sp.expand(p))
print("Accuracy: %e" % np.max([np.abs(f(value) - p.subs('x', value)) for value in xRange]))
plt.plot(x, y, 'ro')
plt.plot(xRange, [p.subs(xSymbol, value) for value in xRange])
plt.show()
