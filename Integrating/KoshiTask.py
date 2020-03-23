import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

eps = 1e-7

def NumericalIntegrate(f, x_0, y_0, b, h, iterations):
    # solve in [x_0, b]
    x_grid = np.arange(x_0, b + eps, h)
    y_i = [y_0] * len(x_grid)
    for _ in range(iterations):
        f_values = [f(x_grid[j], y_i[j]) for j in range(len(x_grid))]

        y_i[1] = y_0 + h / 12 * (5 * f_values[0] + 8 * f_values[1] - f_values[2])
        for j in range(2, len(x_grid)):
            y_i[j] = y_i[j - 2] + h / 3 * \
                     (f_values[j - 2] + 4 * f_values[j - 1] + f_values[j])

    return y_i

x = sp.S('x')

def int_x0_x(f, x_0):
    temp = sp.integrate(f, x)
    return temp - temp.subs(x, x_0)

def AnalyticalIntegrate(f, x_0, y_0, b, h, iterations):
    # solve in [x_0, b]
    x_grid = np.arange(x_0, b + eps, h)
    y_i = y_0
    for i in range(iterations):
        y_i = y_0 + int_x0_x(f(x, y_i), x_0)
    return [y_i.subs(x, X) for X in x_grid]


du = lambda x, u, h, k: h * u + k * np.cos(k * x) - h * np.sin(k * x)
u = lambda x, C, h, k: C * np.exp(h * x) + np.sin(k * x)

# h = -1, k = 1
# u' = f(x, u) = -u + cos(x) + sin(x)
# u(x_0) = y_0
# u_{i} = y_0 + S f(x, u_{i-1}) dx from x_0 to x

h = 0.1
x_grid = np.arange(0, 30 + eps, h)

plt.grid()
plt.plot(x_grid * 1.05 - 1, np.zeros(len(x_grid)), 'black', linewidth=2)
plt.plot([0, 1e-15], [-2, 2], 'black', linewidth=2)
plt.ylim(-2, 2)
plt.xlim(-0.5)

plt.plot(x_grid, AnalyticalIntegrate(lambda x, y: -y + sp.cos(x) + sp.sin(x),
                                     0, sp.S('0'), 30, h, 100),
         'gray', label='analytical integrate', linewidth=5)

plt.plot(x_grid, NumericalIntegrate(lambda x, y: -y + np.cos(x) + np.sin(x),
                                    0, 0, 30, h, 65),
         'lightblue', label='numerical integrate', linestyle='-.')

plt.legend()
plt.show()

