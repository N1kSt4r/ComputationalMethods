import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

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

def f(x):
    return x * sp.tan(x) - 1/3

def fi(x):
    return sp.atan(1 / (3 * x))

def df(x):
    return x*(np.tan(x)**2 + 1) + np.tan(x)

def dfi(x):
    return -1/(3*x**2*(1 + 1/(9*x**2)))

x = sp.Symbol('x')
print('f(x)   = %s' % f(x))
print("f'(x)  = %s" % sp.diff(f(x), x))
print('fi(x)  = tan(0.58*x + 0.1)**0.5')
print("fi'(x) = %s\n" % sp.diff(fi(x), x))

print('Stage of root separation:')
a = 0.1
b = 2
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

# to confirm the monotony
xRange = np.arange(a, b, step=0.0001)
mono = checkMono(df(xRange))
print('Function is monotony on the segment? %s\n' % mono)

if mono:
    flag = True
    x = a
    m = None
    delta = b - a
    q = max(np.abs(dfi(xRange)))

    while flag:
        m = np.abs(fi(x) - x)
        x += eps
        flag = not (m / (1 - q) <= delta)

    if m / (1 - q) <= delta and (0 <= q < 1):
        print('x = fi(x)')
        print("Exist only one root")
        print("Initial approximation is %f" % x)
        print("|x - fi(x)| <= m = %f" % m)
        print("delta = %f" % delta)
        print("q = %f" % q)
        print('m / (1 - q) <= delta ? %s\n' % (m / (1 - q) <= delta))

        k = sp.log(eps / np.abs(m) * (1 - q)) / sp.log(q) // 1 + 1
        print('Theoretical number of steps is %d' % k)

        F = lambda x: (
            (x * fi(fi(x)) - fi(x)**2)
            / (fi(fi(x)) - 2*fi(x) + x)
        )

        k = 0
        while np.abs(x - F(x)) > eps:
            x = F(x)
            if k % 5 == 0:
                print(' New approximation is', x)
            k = k + 1
        print('Practical number of steps is %d\n' % k)

        print('Root is', x)
        print('f(x) = %e' % f(x))


plt.plot(xRange, dfi(xRange))
plt.show()