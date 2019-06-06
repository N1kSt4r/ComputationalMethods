import scipy as sp
from scipy.misc import derivative
import numpy as np
import sympy
import matplotlib.pyplot as plt

# Тема 3, Задание 2, Вариант 21 (7)
# 2.7 x(t) - integrate x(s)ds from 0 to t = t^3/3 - 2t

# x(s)
xSymbol = sympy.Symbol('x')
sSymbol = sympy.Symbol('s')
tSymbol = sympy.Symbol('t')
a = 0
b = 1
down = 0
up = tSymbol
eps = 1e-3

kernel = lambda s, t: 1
def F(x):
    f = sympy.integrate(kernel(sSymbol, tSymbol) * x, sSymbol)
    return (f.subs(sSymbol, up) - f.subs(sSymbol, down)).subs(tSymbol, sSymbol) + sSymbol**3/3 - 2 * sSymbol
norm = lambda x: max(np.abs(x.subs(sSymbol, value)) for value in np.linspace(a, b, 100))

intKernel = sympy.integrate(kernel(sSymbol, tSymbol), sSymbol)
alfa = max(
    np.abs((intKernel.subs(sSymbol, up) - intKernel.subs(sSymbol, down))
        .subs(tSymbol, value)) for value in np.linspace(a + eps, b - eps, 100))

print("Ядро интеграла:", kernel(sSymbol, tSymbol))

x = 0
while max(np.abs((F(x) - x).subs(sSymbol, value)) for value in np.linspace(a, b, 100)) > eps:
    print("", x)
    x = F(x)

print("Решение уравнения: ")
print(" x(t) =", x.subs(sSymbol, tSymbol).evalf(5))
print(" x(t) - Sx(s)ds =",(
    x.subs(sSymbol, tSymbol)
      - sympy.integrate(kernel(sSymbol, tSymbol) * x, sSymbol).subs(sSymbol, up)
      + sympy.integrate(kernel(sSymbol, tSymbol) * x, sSymbol).subs(sSymbol, down)
    ).evalf(5)
)
