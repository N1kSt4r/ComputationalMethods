import numpy as np
from scipy.misc import derivative
from scipy import integrate
from scipy.special import legendre

eps = 1e-6
func = lambda x: np.exp(-x ** 2) * np.sin(x) / (1 + x ** 2)
a, b = [0, 1]

def highest_accuracy_quadratic_formula(func, a, b, n):
    P = legendre(n + 1)
    x = np.roots(P)
    I = 0
    shift_func = lambda x: (b - a)/2 * func(((b - a)/2) * x + (b + a)/2)
    for k in range(n + 1):
        Ak = 2 / (1 - x[k]**2) / derivative(P, x[k], dx=1e-7)**2
        print(f' A[{k}] = {Ak}')
        I += Ak * shift_func(x[k])
    print(' integral by highest accuracy quadratic formula:', I)

print('Highest accuracy quadratic formula')
highest_accuracy_quadratic_formula(func, a, b, 3)
print(f'\nTrue value of integral: {integrate.quad(func, a, b)[0]}')
