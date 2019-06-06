import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def myPrint(obj, title = ""):
    print(title)
    try:
        obj.shape[1]
        for row in obj:
            for elem in row:
                print("%10.7f" % elem, end=' ')
            print()
    except IndexError:
        for elem in obj:
            print("%10.7f" % elem, end=' ')
        print()


def getNextCoef(coef):
    n = len(coef)
    newCoef = np.zeros(n)
    for i in range(n):
        newCoef[i] = coef[i] * coef[i]
        for j in range(1, min(i + 1, n - i)):
            newCoef[i] += 2 * (-1)**j * coef[i + j] * coef[i - j]
        newCoef[i] *= (-1) ** i
    return newCoef

pol = np.poly([1, -5, 9, 14])
a = np.array(pol)

xtemp = np.zeros(len(a) - 1)
x = np.zeros(len(a) - 1)
for i in range(len(x)):
    x[i] = abs(a[i + 1] / a[i])

eps = 1e-5

myPrint(a, "0 coefs:")
myPrint(x, "0 aproximation:")

k = 0
while np.linalg.norm(x - xtemp) > eps:
    xtemp = np.copy(x)

    k += 1
    a = getNextCoef(a)
    a /= np.linalg.norm(a)
    for i in range(len(x)):
        x[i] = (abs(a[i + 1] / a[i])) ** (1 / (2**k))

    myPrint(a, "%d coefs:" % k)
    myPrint(x, "%d aproximation:" % k)

print()
for root in x:
    print(root if abs(np.polyval(pol, root)) < abs(np.polyval(pol, -root)) else -root)

