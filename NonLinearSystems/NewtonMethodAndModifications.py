import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def myFunc(x):
    return sp.tan(0.58*x + 0.1) - x**2

def f(x):
    return np.tan(0.58*x + 0.1) - x**2

def df(x):
    return -2*x + 0.58 * np.tan(0.58*x + 0.1)**2 + 0.58

def d2f(x):
    return 0.58 * (1.16 * np.tan(0.58*x + 0.1)**2 + 1.16) * np.tan(0.58*x + 0.1) - 2

def noZeroValue(values):
    for value in values:
        if np.abs(value) < 1e-15:
            return False
    return True

def checkInitialApproximation(a, b, eps):
    value = a
    h = -f(value) / df(value)
    S = np.arange(min(value, value + 2*h), max(value, value + 2*h), step=eps)
    M = np.max(np.abs(d2f(S)))
    if 2 * np.abs(h) * M <= np.abs(df(value)) and f(value) * d2f(value) >= 0:
        print("h = %f" % h)
        print("max|f''(x)| = %f" % M)
        print("Initial approximation is %f\n" % value)
        q = 0
        for x in S:
            q = max(q, np.abs(d2f(x) / (2 * df(x))))
        k = np.log2(np.log(q * eps) / np.log(q * np.abs(f(value) / df(value))))
        print('Theoretical number of steps is %d' % (int(k) + 1))
        return value
    value = b
    h = -f(value) / df(value)
    S = np.arange(min(value, value + 2*h), max(value, value + 2*h), step=eps)
    M = np.max(np.abs(d2f(S)))
    if 2 * np.abs(h) * M <= np.abs(df(value)) and f(value) * d2f(value) >= 0:
        print("h = %f" % h)
        print("max|f''(x)| = %f" % M)
        print("Initial approximation is %f\n" % value)
        q = 0
        for x in S:
            q = max(q, np.abs(d2f(x) / (2 * df(x))))
        k = np.log2(np.log(q * eps) / np.log(q * np.abs(f(value) / df(value))))
        print('Theoretical number of steps is %d' % (int(k) + 1))
        return value

x = sp.Symbol('x')
print('f(x)   = %s' % myFunc(x))
print("df(x)  = %s" % sp.diff(myFunc(x), x))
print("d2f(x) = %s" % sp.diff(sp.diff(myFunc(x), x),x))

xRange = np.arange(0, 2, step=0.01)
plt.plot(xRange, 0 * xRange, 'g--')
plt.plot(xRange, f(xRange), 'red')
plt.plot(xRange, df(xRange), 'orange')
plt.plot(xRange, d2f(xRange), 'yellow')
plt.show()

print('\nStage of root separation:')
a = 0.1
b = 100
eps = 1e-5

while b - a > 0.1:
    print(' a = %f, f(a) = %f' % (a, f(a)))
    print(' b = %f, f(b) = %f' % (b, f(b)))
    middle = (b + a) / 2
    print(' length = %f, f(middle) = %f\n' % (b - a, f(middle)))
    if f(a) * f(middle) < 0:
        b = middle
    else:
        a = middle

print('Final ends of segment are a = %.4f, b = %.4f\n' % (a, b))

xRange = np.arange(a, b, step=eps)

if noZeroValue(df(xRange)) and noZeroValue(d2f(xRange)):
    print("Function's derivative and second derivative haven't zero values")
    print("Exist only one root\n")
    print('x = fi(x)')
    print("fi(x) = x - f(x)/f'(x)")
    x0 = checkInitialApproximation(a, b, eps)

    print("\nNewton's method:")
    xtemp = x0
    x = xtemp - f(xtemp)/df(xtemp)
    k = 0
    while np.abs(x - xtemp) > eps:
        k += 1
        print(" %d approximation - %f" % (k, x))
        xtemp = x
        x = xtemp - f(xtemp)/df(xtemp)

    print('Practical number of steps is %d' % k)
    print('Root is', x)
    print('f(x) = %e' % f(x))

    print("\nNewton's method with const derivative:")
    xtemp = x0
    dfx0 = df(x0)
    x = xtemp - f(xtemp)/dfx0;
    k = 0
    while np.abs(x - xtemp) > eps:
        k += 1
        print(" %d approximation - %f" % (k, x))
        xtemp = x
        x = xtemp - f(xtemp)/dfx0

    print('Practical number of steps is %d' % k)
    print(' Root is', x)
    print(' f(x) = %e' % f(x))

    print("\nMethod of secants:")
    xtemptemp = x0 - f(x0)/df(x0)
    xtemp = xtemptemp - f(xtemptemp)/df(xtemptemp)
    print(" %d approximation - %f" % (1, xtemptemp))
    print(" %d approximation - %f" % (2, xtemp))
    dfx0 = df(x0)
    x = xtemp - f(xtemp) * (xtemp - xtemptemp) / (f(xtemp) - f(xtemptemp))
    k = 2
    while np.abs(x - xtemp) > eps:
        k += 1
        print(" %d approximation - %f" % (k, x))
        xtemp = x
        x = xtemp - f(xtemp) * (xtemp - x0) / (f(xtemp) - f(x0))

    print('Practical number of steps is %d' % k)
    print(' Root is', x)
    print(' f(x) = %e' % f(x))

    print("\nMethod of chords:")
    xtemp = x0
    x = a if x0 == b else b
    print(" x0 = %f" % xtemp)
    print(" x1 = %f" % x)
    dfx0 = df(x0)
    k = 0
    while np.abs(x - xtemp) > eps:
        k += 1
        print(" %d approximation - %f" % (k, x))
        xtemp = x
        x = xtemp - f(xtemp) * (xtemp - x0) / (f(xtemp) - f(x0))

    print('Practical number of steps is %d' % k)
    print(' Root is', x)
    print(' f(x) = %e' % f(x))

else:
    print("Function's derivative and second derivative have zero values")
