import scipy as sp
from scipy.misc import derivative
import numpy as np
import sympy
import matplotlib.pyplot as plt

# Тема 2, Задание 2, Вариант 21 (7)
# 2.7 3x^2 - 10x - 141 = 0

xSymbol = sympy.Symbol('x')
eps = 1e-4
g = lambda x: 3*x**2 - 10*x - 141
print("g(x) =", g(xSymbol))
xRange = np.arange(-10, 10, 0.1)
plt.plot(xRange, g(xRange))
plt.plot(xRange, 0 * xRange)
plt.show()


def calculate(g, f):
    absdf = lambda x: np.abs(derivative(f, x, dx=eps))
    print("_____________________________\n")
    print("f(x) =", f(xSymbol))
    print("f'(x) =", sympy.diff(f(xSymbol), xSymbol))

    xRange = np.arange(-10, 10.1, 0.1)
    xRange = xRange[np.where(np.abs(absdf(xRange)) < 1)]
    if len(xRange) > 1:
        print("Отображение сжимающее для x из (%.4f, %.4f)" % (xRange[0], xRange[len(xRange) - 1]))

    found = False
    for a in xRange:
        for r0 in np.arange(0.1, 0.5, 0.1):
            L = max(absdf(np.arange(a - r0, a + r0, eps)))
            if np.abs(a - f(a)) < (1 - L) * r0 and L < 1:
                print("Сжимающий коэффициент: %.4f" % L)
                print("%.4f < %.4f" % (np.abs(a - f(a)), (1 - L) * r0))
                print("%.4f < 1" % L)
                print("Условия выполняются")
                print("a = %.4f, r0 = %.4f" % (a, r0))
                estimate = int(np.ceil(np.log((1 - L) * eps / np.abs(a - f(a))) / np.log(L)))
                found = True
                print("Примерная оценка: %d" % estimate)
                for i in range(estimate):
                    a = f(a)
                print("Значение g(%.5f) = %.4e" % (a, g(a)))
                break
        if found:
            break
    if not found:
        print("Условия не выполняются")


f = lambda x: 0.3*x**2 - 14.1
calculate(g, f)

f = lambda x: ((10*x + 141) / 3)**0.5
calculate(g, f)

f = lambda x: 141 / (3*x - 10)
calculate(g, f)
