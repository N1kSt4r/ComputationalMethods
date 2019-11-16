import numpy as np
from scipy.misc import derivative
from scipy import integrate

eps = 1e-6
func = lambda x: np.exp(-x ** 2) * np.sin(x) / (1 + x ** 2)
a, b = [0, 1]


def next_Runge(I1, I2, k):
    return (I2 - I1) / (1 - 2**-k)


def Simpson_integrate(func, a, b, N):
    h = (b - a) / N
    func_values = func(np.linspace(a, b, N + 1))
    I = func_values[0] + func_values[N]
    I += np.sum(func_values[2:N:2]) * 2
    I += np.sum(func_values[1:N:2]) * 4
    return I * h / 3


def Simpson_Runge(func, a, b):
    max_d4f = max([abs(derivative(func, x, dx=0.0001, order=5, n=4)) for x in np.linspace(a, b, 100)])
    print(f'  max(abs(d4f)) on [{a}, {b}]: {max_d4f}')
    N = int(((b - a) ** 5 * max_d4f / 180 / eps) ** (1 / 4)) // 2 * 2 + 2
    print(f'  a priory estimate of segments: {N}')
    I1 = Simpson_integrate(func, a, b, 10)
    I2 = Simpson_integrate(func, a, b, 20)
    print('  integral with step {:.2e}: {}'.format((b - a) / 10, I1))
    print('  integral with step {:.2e}: {}'.format((b - a) / 20, I2))

    N = 20
    while abs(next_Runge(I1, I2, 4)) > eps:
        N *= 2
        I1, I2 = I2, Simpson_integrate(func, a, b, N)
        print('  integral with step {:.2e}: {}'.format((b - a) / N, I2))
    print('  itegral with Runge:', I1 + next_Runge(I1, I2, 4))
    print('  true number of segments:', N)


def trapeze_integrate(func, a, b, N):
    h = (b - a) / N
    func_values = func(np.linspace(a, b, N + 1))
    I = (func_values[0] + func_values[N]) / 2
    I += np.sum(func_values[1:-1])
    return I * h


def trapeze_Runge(func, a, b):
    max_d2f = max([abs(derivative(func, x, dx=0.0001, order=5, n=2)) for x in np.linspace(a, b, 100)])
    print(f'  max(abs(d2f)) on [{a}, {b}]: {max_d2f}')
    N = int(((b - a) ** 3 * max_d2f / 12 / eps) ** (1 / 2)) // 2 * 2 + 2
    print(f'  a priory estimate of segments: {N}')
    I1 = trapeze_integrate(func, a, b, 10)
    I2 = trapeze_integrate(func, a, b, 20)
    print('  integral with step {:.2e}: {}'.format((b - a) / 10, I1))
    print('  integral with step {:.2e}: {}'.format((b - a) / 20, I2))

    N = 20
    while abs(next_Runge(I1, I2, 2)) > eps:
        N *= 2
        I1, I2 = I2, trapeze_integrate(func, a, b, N)
        print('  integral with step {:.2e}: {}'.format((b - a) / N, I2))
    print('  next approximation with Runge: {}'.format(I1 + next_Runge(I1, I2, 2)))
    print('  true number of segments:', N)


print('Simpson')
Simpson_Runge(func, a, b)
print('\nTrapeze')
trapeze_Runge(func, a, b)
print(f'\nTrue value of integral: {integrate.quad(func, a, b)[0]}')
