import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

a, b = -1, 4

values = [-1, 0, 2, 3]
fvalue = [0, 2, 1, 3]

basis = [(lambda power: lambda x: x**power)(power) for power in range(6)]

A = np.array([
    [
        sum([x(value) * y(value) for value in values]) for x in basis
    ] for y in basis
])
B = np.array([
    sum([x(value) * fvalue[i] for i, value in enumerate(values)]) for x in basis
])
coefs = np.linalg.solve(A, B)
print(A, B, coefs)

def g(x):
    sum = 0
    for i in range(len(coefs)):
        sum += coefs[i] * basis[i](x)
    return sum

xRange = np.linspace(a, b, (b - a) * 100)
plt.plot(xRange, g(xRange))
plt.show()