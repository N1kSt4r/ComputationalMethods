import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def myPrint(obj, title="", format="f"):
    print(title)
    if len(obj.shape) > 1:
        for row in obj:
            for elem in row:
                print(("%10.4" + format) % elem, end=' ')
            print()
    else:
        for elem in obj:
            print(("%10.4" + format) % elem, end=' ')
        print()

eps = 1e-5

f = [
    lambda x: x[0] + x[0]**2 - 2*x[1]*x[2] - 0.1,
    lambda x: x[1] - x[1]**2 + 3*x[0]*x[2] + 0.2,
    lambda x: x[2] + x[2]**2 + 2*x[0]*x[1] - 0.3,
]

def F(x):
    res = np.zeros(len(f))
    for i in range(len(f)):
        res[i] = f[i](x)
    return res

x = [sp.Symbol('x[0]'), sp.Symbol('x[1]'), sp.Symbol('x[2]')]
for i in range(len(x)):
    for j in range(len(x)):
        print("%13s" % sp.diff(f[i](x), x[j]), end='')
    print()

def getJacobi(x):
    return np.array([
        [1 + 2*x[0],    -2*x[2],    -2*x[1]],
        [    3*x[2], 1 - 2*x[1],     3*x[0]],
        [    2*x[1],     2*x[0], 1 + 2*x[2]]
    ])

def nextAprox(x):
    J = getJacobi(x)
    JJf = np.dot(np.dot(J, J.T), F(x))
    t = np.dot(F(x), JJf) / np.dot(JJf, JJf)
    return x - t * np.dot(J.T, F(x))

x = np.array([0, 0, 0])
xtemp = np.copy(x)
x = nextAprox(x)
myPrint(x, "\n1 approximation:")
k = 2

while np.max(np.abs(x - xtemp)) > eps:
    xtemp = np.copy(x)
    x = nextAprox(x)
    myPrint(x, "%d approximation:" % k)
    k += 1

myPrint(F(x), "F(x):", "e")
