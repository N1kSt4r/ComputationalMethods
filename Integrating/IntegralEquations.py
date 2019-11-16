import numpy as np


def display(result):
    for num in result:
        print('%10.5f' % num, end=' ')
    print()


def Simpson_integrate(func_values, step):
    I = func_values[0] + func_values[-1]
    I += np.sum(func_values[2:-1:2]) * 2
    I += np.sum(func_values[1:-1:2]) * 4
    return I * step / 3


def mechanic_quad_fredholm(K, f, h, a, b, n):
    points = np.linspace(a, b, n + 1)
    step = (b - a) / n
    A = np.ones(n + 1) * step
    A[[0, n]] /= 2

    M = np.identity(n + 1) + -h * A * np.array([[K(x, s) for s in points] for x in points])
    return np.linalg.solve(M, f(points))


def sequentional_approximations_fredholm(K, f, h, a, b, n=10):
    n *= 100
    x = np.linspace(a, b, n + 1)
    K_values = np.array([K(x, s) for s in x])
    u_prev = np.ones(n + 1)
    u = np.zeros(n + 1)

    while np.linalg.norm(u - u_prev) > 1e-5:
        u_prev = u.copy()
        for i in range(n + 1):
            u[i] = f(x[i]) + h * Simpson_integrate(u * K_values[i], (b - a) / n)

    return u[::100]


def get_alfa(points, index, value):
    numer = np.prod(value - np.append(points[:index], points[index + 1:]))
    denum = np.prod(points[index] - np.append(points[:index], points[index + 1:]))
    return numer / denum


def kernel_replacement_by_degenerate(K, f, h, a, b, n):
    pivots = np.linspace(a, b, n + 1)
    values = np.linspace(a, b, 100 * n + 1)

    alfa_values = np.array([[get_alfa(pivots, i, value)
                             for value in values] for i in range(n + 1)])
    beta_values = np.array([[K(pivots[i], value)
                             for value in values] for i in range(n + 1)])

    a_values = np.array([[Simpson_integrate(beta_values[i] * alfa_values[j], (b - a) / n / 100)
                          for j in range(n + 1)] for i in range(n + 1)])
    b_values = np.array([Simpson_integrate(beta_values[i] * f(values), (b - a) / n / 100)
                         for i in range(n + 1)])

    M = np.identity(n + 1) - h * a_values
    c = np.linalg.solve(M, b_values)
    return f(pivots) + h * c


def mechanic_quad_volter(K, f, h, a, b, n):
    x = np.linspace(a, b, n + 1)
    step = (b - a) / n
    A = np.ones(n + 1) * step
    A[[0, n]] /= 2

    M = np.array([[0 if i < j else K(x[i], x[j]) for j in range(n + 1)] for i in range(n + 1)])
    M = np.identity(n + 1) + -h * A * M
    u = np.linalg.solve(M, f(x))
    return u


# 1. u(x) - 0.5 S exp(-xy) u(y) dy from 0 to 1 = x
K = lambda x, s: np.exp(-x * s)
F = lambda x: x
h = 0.5
a = 0
b = 1
n = 10

print('Mechanic quad Fredholm:')
display(mechanic_quad_fredholm(K, F, h, a, b, n))
print('Sequentional approximations Fredholm:')
display(sequentional_approximations_fredholm(K, F, h, a, b, n))
print('Kernel replacement by degenerate')
display(kernel_replacement_by_degenerate(K, F, h, a, b, n))

# 2. u(x) - S sin(x - y) u(y) dy from 0 to x = cos(x)  x from 0 to pi/2
K = lambda x, s: np.sin(x - s)
F = lambda x: np.cos(x)
h = 1
a = 0
b = np.pi / 2
n = 10

print('Mechanic quad Volter:')
display(mechanic_quad_volter(K, F, h, a, b, n))
