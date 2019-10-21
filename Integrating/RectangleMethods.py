import numpy as np
from scipy.misc import derivative
from scipy import integrate

eps = 1e-6
func = lambda x: np.exp(-x ** 2) * np.sin(x) / (1 + x ** 2)
a, b = [0, 1]


def next_Runge(I1, I2, k):
    return (I2 - I1) / (1 - 2**-k)


def left_rectangle_integrate(func, a, b, N):
    h = (b - a) / N
    func_values = list(map(func, np.linspace(a, b - h, N)))
    I = sum(func_values)
    return I * h


def left_rectangle_Runge(func, a, b):
    max_df = max([abs(derivative(func, x, dx=1e-6)) for x in np.linspace(a, b, 100)])
    print(f'  max(abs(df)) on [{a}, {b}]: {max_df}')
    N = int((b - a) ** 2 * max_df / 2 / eps) // 2 * 2 + 2
    print(f'  a priory estimate of segments: {N}')
    I1 = left_rectangle_integrate(func, a, b, 10)
    I2 = left_rectangle_integrate(func, a, b, 20)
    print('  integral with step {:.2e}: {}'.format((b - a) / 10, I1))
    print('  integral with step {:.2e}: {}'.format((b - a) / 20, I2))

    N = 20
    true_iterations = 2
    while abs(next_Runge(I1, I2, 1)) > eps:
        N *= 2
        true_iterations += 1
        I1, I2 = I2, left_rectangle_integrate(func, a, b, N)
        print('  integral with step {:.2e}: {}'.format((b - a) / N, I2))
    print('  next approximation with Runge: {}'.format(I1 + next_Runge(I1, I2, 2)))
    print('  true number of segments:', 10 * 2**(true_iterations - 1))


def right_rectangle_integrate(func, a, b, N):
    h = (b - a) / N
    func_values = list(map(func, np.linspace(a + h, b, N)))
    I = sum(func_values)
    return I * h


def right_rectangle_Runge(func, a, b):
    max_df = max([abs(derivative(func, x, dx=1e-6)) for x in np.linspace(a, b, 100)])
    print(f'  max(abs(df)) on [{a}, {b}]: {max_df}')
    N = int((b - a) ** 2 * max_df / 2 / eps) // 2 * 2 + 2
    print(f'  a priory estimate of segments: {N}')
    I1 = right_rectangle_integrate(func, a, b, 10)
    I2 = right_rectangle_integrate(func, a, b, 20)
    print('  integral with step {:.2e}: {}'.format((b - a) / 10, I1))
    print('  integral with step {:.2e}: {}'.format((b - a) / 20, I2))

    N = 20
    true_iterations = 2
    while abs(next_Runge(I1, I2, 1)) > eps:
        N *= 2
        true_iterations += 1
        I1, I2 = I2, right_rectangle_integrate(func, a, b, N)
        print('  integral with step {:.2e}: {}'.format((b - a) / N, I2))
    print('  next approximation with Runge: {}'.format(I1 + next_Runge(I1, I2, 2)))
    print('  true number of segments:', 10 * 2**(true_iterations - 1))


def middle_rectangle_integrate(func, a, b, N):
    h = (b - a) / N
    func_values = list(map(func, np.linspace(a + h / 2, b - h / 2, N)))
    I = sum(func_values[1:])
    return I * h


def middle_rectangle_Runge(func, a, b):
    max_df = max([abs(derivative(func, x, dx=1e-6)) for x in np.linspace(a, b, 100)])
    print(f'  max(abs(df)) on [{a}, {b}]: {max_df}')
    N = int((b - a) ** 2 * max_df / 2 / eps) // 2 * 2 + 2
    print(f'  a priory estimate of segments: {N}')
    I1 = middle_rectangle_integrate(func, a, b, 10)
    I2 = middle_rectangle_integrate(func, a, b, 20)
    print('  integral with step {:.2e}: {}'.format((b - a) / 10, I1))
    print('  integral with step {:.2e}: {}'.format((b - a) / 20, I2))

    N = 20
    true_iterations = 2
    while abs(next_Runge(I1, I2, 1)) > eps:
        N *= 2
        true_iterations += 1
        I1, I2 = I2, middle_rectangle_integrate(func, a, b, N)
        print('  integral with step {:.2e}: {}'.format((b - a) / N, I2))
    print('  next approximation with Runge: {}'.format(I1 + next_Runge(I1, I2, 2)))
    print('  true number of segments:', 10 * 2**(true_iterations - 1))


print('\nLeft rectangle')
left_rectangle_Runge(func, a, b)
print('\nRight rectangle')
right_rectangle_Runge(func, a, b)
print('\nMiddle rectangle')
middle_rectangle_Runge(func, a, b)
print(f'\nTrue value of integral: {integrate.quad(func, a, b)[0]}')

