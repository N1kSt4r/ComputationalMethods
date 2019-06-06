import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def ChebyshevNodes(a, b, n):
    x = [(a + b) / 2] * n
    for k in range(n):
        x[k] -= (b - a)*np.cos((k + 0.5) / n * np.pi) / 2
    return x

n = 5
a, b = 0, np.pi / 2
xSymbol = sp.S('x')
f = lambda x: np.sin(x + 0.5)
x = np.array(ChebyshevNodes(a, b, n))
print(x)
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

xRange = np.linspace(a, b, 100)
print("P(x) =", sp.expand(p))
print("Accuracy: %e" % np.max([np.abs(f(value) - p.subs('x', value)) for value in xRange]))
plt.plot(x, y, 'ro')
plt.plot(xRange, [p.subs(xSymbol, value) for value in xRange])
plt.show()
