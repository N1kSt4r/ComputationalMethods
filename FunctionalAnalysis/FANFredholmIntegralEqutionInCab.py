import scipy as sp
from scipy.misc import derivative
import numpy as np
import sympy
import matplotlib.pyplot as plt

# Тема 3, Задание 1, Вариант 21 (7)
# 1.21 a = 0; b = 1; x(t) - h * integrate s^2(t^2 + 1)x(s) ds from 0 to 1 = 1

# x(s)
xSymbol = sympy.Symbol('x')
sSymbol = sympy.Symbol('s')
tSymbol = sympy.Symbol('t')
a = 0
b = 1
down = 0
up = 1
eps = 1e-3
h = 1.11

kernel = lambda s, t: s**2 * (t**2 + 1)
def F(x):
    f = sympy.integrate(kernel(sSymbol, tSymbol) * x, sSymbol)
    return h * (f.subs(sSymbol, up) - f.subs(sSymbol, down)).subs(tSymbol, sSymbol) + 1
norm = lambda x: np.abs(h) * max(np.abs(x.subs(sSymbol, value)) for value in np.linspace(a, b, 100))

intKernel = sympy.integrate(kernel(sSymbol, tSymbol), sSymbol)
alfa = max(
    np.abs((intKernel.subs(sSymbol, up) - intKernel.subs(sSymbol, down))
        .subs(tSymbol, value)) for value in np.linspace(a, b, 100))

print("Ядро интеграла:", kernel(sSymbol, tSymbol))
print("Коэффициент сжатия: %.3f" % alfa)
print(" => для |h| < %.3f можно применить метод сжимающих отображений" % (1/alfa))

x = sympy.S(0)
print("||x0 - x1|| = %.2f" % norm(x - F(x)))
estimate = int(np.log(float(eps * (1 - alfa) / norm(x - F(x)))) / np.log(float(alfa))) + 1
print("Примерная оценка: %d" % estimate)
while norm(x - F(x)) > eps:
    print("", x.evalf(5))
    x = F(x)

print("Решение уравнения при h = %.2f" % h)
print(" x(t) =", x.subs(sSymbol, tSymbol).evalf(5))
print(" x(t) - %.2f Ss^2(t^2 + 1)x(s)ds =" % h,(
    x.subs(sSymbol, tSymbol)
      - h * sympy.integrate(kernel(sSymbol, tSymbol) * x, sSymbol).subs(sSymbol, up)
      + h * sympy.integrate(kernel(sSymbol, tSymbol) * x, sSymbol).subs(sSymbol, down)
    ).evalf(5)
)
