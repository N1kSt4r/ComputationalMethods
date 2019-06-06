import sympy as sp
import numpy as np

f = [
    lambda x: sp.sin((x[0] + x[1]) / 3) - 2*x[0],
    lambda x: sp.cos((x[0] - x[1]) / 3) - 2*x[1]
]

def getJacobi(f, values):
    Jacobi = np.zeros((len(values), len(values)))
    variables = [sp.Symbol("x[%d]" % index) for index in range(len(values))]

    for i in range(len(values)):
        for j in range(len(values)):
            Jacobi[i][j] = sp.diff(f[i](variables), variables[j])\
                .subs(list(zip(variables, values)))
    return Jacobi

eps = 1e-3
epsInner = 1e-5
x = np.array([0.1, 0.4])

nextApproximationInner = lambda x, i: x[i] - f[i](x) / getJacobi(f, x)[i][i]

def nextApproximationOuter(vector):
    x = np.copy(vector)
    for i in range(len(f)):
        while np.abs(x[i] - nextApproximationInner(x, i)) > epsInner:
            x[i] = nextApproximationInner(x, i)
    return x

k = 0
while np.max(np.abs(x - nextApproximationOuter(x))) > eps:
    print("x(%d) = " % k, ("%11.6f" * len(x)) % tuple(value for value in x))
    x = nextApproximationOuter(x)
    k += 1

print("x    = ", ("%11.6f" * len(x)) % tuple(value for value in x))
print("F(x) = ", ("%.4e  " * len(f)) % tuple(f[i](x) for i in range(len(f))))
