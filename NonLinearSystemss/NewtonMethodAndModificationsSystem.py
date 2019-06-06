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

def F(x):
    res = np.zeros(len(f))
    for i in range(len(f)):
        res[i] = f[i](x)
    return res

def getJacobi(f, values):
    Jacobi = np.zeros((len(values), len(values)))
    variables = [sp.Symbol("x[%d]" % index) for index in range(len(values))]

    for i in range(len(values)):
        for j in range(len(values)):
            Jacobi[i][j] = sp.diff(f[i](variables), variables[j])\
                .subs(list(zip(variables, values)))
    return Jacobi

eps = 1e-3

print("\nNewton method:")
x = np.array([0.1, 0.4])
xtemp = np.copy(x)
x = x - np.dot(np.linalg.inv(getJacobi(f, x)), F(x))
myPrint(x, "1 approximation:")
k = 2

while np.max(np.abs(x - xtemp)) > eps:
    xtemp = np.copy(x)
    x = x - np.dot(np.linalg.inv(getJacobi(f, x)), F(x))
    myPrint(x, "%d approximation:" % k)
    k += 1

myPrint(F(x), "F(x):", "e")

print("\nNewton method with const Jacobi:")
x = np.array([0.1, 0.4])
xtemp = np.copy(x)
x = x - np.dot(np.linalg.inv(getJacobi(f, x)), F(x))
myPrint(x, "1 approximation:")
k = 2

J = np.linalg.inv(getJacobi(f, x))
while np.max(np.abs(x - xtemp)) > eps:
    xtemp = np.copy(x)
    x = x - np.dot(J, F(x))
    myPrint(x, "%d approximation:" % k)
    k += 1

myPrint(F(x), "F(x):", "e")

def getEqualJacobi(x, xtemp):
    n = len(x)
    J = np.zeros(n**2).reshape(n, n)
    for i in range(n):
        for j in range(n):
            xDiffJ = np.copy(x)
            xDiffJ[j] = xtemp[j]
            J[i][j] = (f[i](x) - f[i](xDiffJ)) / (x[j] - xtemp[j])
    return J

print("\nSecants method:")
x = np.array([0.1, 0.4])
xtemp = np.copy(x)
x = x - np.dot(np.linalg.inv(getJacobi(f, x)), F(x))
myPrint(x, "1 approximation:")
k = 2

while np.max(np.abs(x - xtemp)) > eps:
    xtemp, x = np.copy(x), x - np.dot(np.linalg.inv(getEqualJacobi(x, xtemp)), F(x))
    myPrint(x, "%d approximation:" % k)
    k += 1

myPrint(F(x), "F(x):", "e")
