import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def myPrint(obj, title="", format="f"):
    print(title)
    if len(obj.shape) > 1:
        for row in obj:
            for elem in row:
                print(("%10.6" + format) % elem, end=' ')
            print()
    else:
        for elem in obj:
            print(("%10.6" + format) % elem, end=' ')
        print()

f = [
    lambda x: sp.sin((x[0] + x[1]) / 3) - 2*x[0],
    lambda x: sp.cos((x[0] - x[1]) / 3) - 2*x[1]
]

fi = [
    lambda x: sp.sin((x[0] + x[1]) / 3) / 2,
    lambda x: sp.cos((x[0] - x[1]) / 3) / 2
]

print("\nSimple iteration method:")
x = np.array([0.1, 0.4], float)
myPrint(x, "initial approximation:")

eps = 1e-3
k = 0

def nextApproximation(x):
    next = np.copy(x)
    for i in range(len(x)):
        next[i] = fi[i](x)
    return next

while np.linalg.norm(x - nextApproximation(x)) > eps:
    x = nextApproximation(x)
    k += 1
    myPrint(x, "%d approximation: " % k)

print("\nF(x): %10.4e %10.4e" % (f[0](x), f[1](x)))


print("\nZeidel method:")
x = np.array([0.1, 0.4])
myPrint(x, "initial approximation:")

xtemp = np.copy(x)
for i in range(len(x)):
    x[i] = fi[i](xtemp)

eps = 1e-3
k = 0

while np.linalg.norm(x - xtemp) > eps:
    k += 1
    xtemp = np.copy(x)
    for i in range(len(x)):
        x[i] = fi[i](x)
    myPrint(x, "%d approximation: " % k)

print("\nF(x): %10.4e %10.4e" % (f[0](x), f[1](x)))