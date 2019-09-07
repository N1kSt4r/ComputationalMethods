import numpy as np
import scipy as sp
from scipy.misc import derivative
from scipy import integrate
import matplotlib.pyplot as plt

w = lambda arg: np.prod(arg - x)
dw = lambda arg, k: np.prod(arg - np.append(x[:k], x[k+1:]))

a, b = 1, 2
f = lambda x: 0.8 * np.log(x) + 0.6 * np.cos(x)
max_df = max([abs(derivative(f, x, dx=0.0001)) for x in np.linspace(a, b, 100)])
max_d2f = max([abs(derivative(f, x, dx=0.0001, n=2)) for x in np.linspace(a, b, 100)])
eps = 1e-7

print("Scipy value:           ", integrate.quad(f, a, b)[0])

# left rectangle method

n = int(max_df * (b - a)**2 / 4/eps) + 1
dx = (b - a) / n
x = a
value = 0

for i in range(n):
    value += dx * f(x)
    x += dx
print("Left rectangle method: ", value)

# right rectangle method

n = int(max_df * (b - a)**2 / 4/eps) + 1
dx = (b - a) / n
x = b
value = 0

for i in range(n):
    value += dx * f(x)
    x -= dx
print("Right rectangle method:", value)

# trapeze method

n = int((max_d2f * (b - a)**3 / 12 / eps)**0.5) + 1
dx = (b - a) / n
x = a + dx/2
value = 0

for i in range(n):
    value += dx * f(x)
    x += dx
print("trapeze method:        ", value)
