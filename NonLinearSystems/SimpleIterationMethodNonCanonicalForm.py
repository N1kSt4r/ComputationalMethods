import sympy as sp
import numpy as np

def checkMono(values):
    if values[0] == 0:
        return False
    if values[0] > 0:
        for value in values:
            if value < 0:
                return False
    else:
        for value in values:
            if value > 0:
                return False
    return True

def less(values, x):
    for value in values:
        if value >= x:
            return False
    return True

def f(x):
    return sp.tan(0.58*x + 0.1) - x**2

def df(x):
    return -2*x + 0.58 * np.tan(0.58*x + 0.1)**2 + 0.58


x = sp.Symbol('x')
print('f(x)   = tan(0.58x + 0.1) - x^2')
print("f'(x)  = %s\n" % sp.diff(f(x), x))

print('Stage of root separation:')
a = 0.1
b = 10
eps = 1e-15

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

xRange = np.arange(a, b, step=0.0001)
mono = checkMono(df(xRange))
print('Function is monotony on the segment? %s\n' % mono)

if mono:
    print('fi(x) = x - f(x)/M')
    print("Exist only one root")
    x = (a + b) / 2
    print("Initial approximation is %f" % x)
    m = np.min(np.abs(df(xRange)))
    M = np.max(np.abs(df(xRange)))

    fi = lambda x: x - f(x) / M
    if df(x) < 0:
        fi = lambda x: x + f(x) / M
    q = 1 - m / M

    print('m = %f, M = %f, q = 1 - m / M = %f\n' % (m, M, q))

    k = sp.log(eps / np.abs(m) * (1 - q)) / sp.log(q) // 1 + 1
    print('Theoretical number of steps is %d' % k)

    k = 0
    while np.abs(x - fi(x)) > eps:
        x = fi(x)
        if k % 5 == 0:
            print(' New approximation is', x)
        k = k + 1
    print('Practical number of steps is %d\n' % k)

    print('Root is', x)
    print('f(x) = %e' % f(x))

