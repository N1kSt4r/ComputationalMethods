import scipy as sp
from scipy.integrate import quad as integrate
import numpy as np
import matplotlib.pyplot as plt

a, b = 1, 2
p = lambda x: 1
f = lambda x: 0.8 * sp.log(x) + 0.6 * sp.cos(x)

basis = [(lambda power: lambda x: x**power)(power) for power in range(6)]

A = np.array([[integrate(lambda value: p(value) * x(value) * y(value), a, b)[0] for x in basis] for y in basis])
B = np.array([integrate(lambda value: p(value) * x(value) * f(value), a, b)[0] for x in basis])
coefs = np.linalg.solve(A, B)

def g(x):
    sum = 0
    for i in range(len(coefs)):
        sum += coefs[i] * basis[i](x)
    return sum

xRange = np.linspace(a, b, (b - a) * 100)
print("Accuracy: %e" % np.max(np.abs(f(xRange) - g(xRange))))
plt.plot(xRange, f(xRange))
plt.plot(xRange, g(xRange))
plt.show()
