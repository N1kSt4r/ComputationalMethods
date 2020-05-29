import numpy as np
from numpy import power, sin, cos, tan
from scipy.linalg import norm, solve_banded

true_ans = np.linspace(-1, 0, 101)
f = lambda x: x * sin(2 * x)
k = lambda x: power(cos(x), 2)
q = lambda x: sin(2 * x)
dk = lambda x: -sin(2 * x)
a_0 = 0
b_0 = -1
a_1 = 1
b_1 = cos(1)**2

x = np.linspace(0, 1, 101)[1: -1]
i_range = np.arange(1, x.size - 1)
h = 0.01
f_v = f(x)
k_v = k(x)
q_v = q(x)
dk_v = dk(x)

print('True solution:')
for i in true_ans:
    print('{:+.5f}'.format(i), end=' ')
print('\n')

def Solve(A_0_0, A_0_1, A_n_nm1, A_n_n,
          A_below, A_in, A_above,
          b_beg, b_mid, b_end):
    A = np.array([
        [0, A_0_1, *A_above],
        [A_0_0, *A_in, A_n_n],
        [*A_below, A_n_nm1, 0],
    ])

    b = [b_beg, *b_mid, b_end]

    solution = solve_banded((1, 1), A, b)
    print('Solution:')
    for i in solution:
        print('{:.5f}'.format(i), end=' ')
    disc = solution - true_ans
    print('\nDiscrepancy:')
    for i in disc:
        print('{:+.1e}'.format(i), end=' ')
    print('\nRelative error: {:.2e}'.format(norm(disc) / norm(solution)))
    print()

def IncreasingOrderMethod():
    print('Increasing order method: ')
    Solve(-k(0) ** 2 / h - dk(0) * a_0 - q(0) * k(0) - a_0 * k(0), k(0) ** 2 / h,
          -k(1) ** 2 / h, k(1) ** 2 / h - dk(1) * a_1 + q(1) * k(1) + a_1 * k(1),
          (2 * k_v - dk_v * h) / h ** 2 / 2,
          -(2 * k_v + q_v * h ** 2) / h ** 2,
          (2 * k_v + dk_v * h) / h ** 2 / 2,
          -f(0) * k(0) + b_0 * (dk(0) - k(0)), -f_v, f(1) * k(1) + b_1 * (dk(1) + k(1)))

def BalanceMethod():
    print('Balance method: ')
    c = np.append(k(x - h / 2), k(1 - h / 2))
    Solve(-(c[0] / h + a_0 + q(h / 4) * h / 2), c[0] / h,
          c[-1] / h, -(c[-1] / h + a_1 + q(1 - h / 4) * h / 2),
          c[:-1], -(q_v * h ** 2 + c[1:] + c[:-1]), c[1:],
          -b_0 - h / 2 * f(h / 4),
          -f_v * h**2,
          -b_1 - h / 2 * f(1 - h / 4))

def TrapezeInt(F, a, b):
    return (b - a) * (F(a) + F(b)) / 2

def RitzMethod():
    print('Ritz method: ')
    x_ext = np.linspace(0, 1, 101)
    c = (k(x_ext[:-1]) + k(x_ext[1:])) / 2
    Solve(-(a_0 + c[0] / h - h * q(0) / 2), c[0] / h,
          c[-1] / h, -(a_1 + c[-1] / h + h * q(1) / 2),
          c[:-1] / h**2, -(c[:-1] + c[1:]) / h**2 - q_v, c[1:] / h**2,
          -b_0 + h * f(0) / 2, -f_v, -b_1 - h * f(1) / 2)

IncreasingOrderMethod()
BalanceMethod()
RitzMethod()

