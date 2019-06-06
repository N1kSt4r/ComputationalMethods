import numpy as np
import matplotlib.pyplot as plt

def binSearch(arr, value):
    l, r = 0, len(arr)
    while (r - l > 1):
        mid = (r + l) // 2
        if value < arr[mid]:
            r = mid
        else:
            l = mid
    return r

# natural spline

x = np.array([1, 2, 3, 4, 5, 6])
y = np.array([1.0002, 1.0341, 0.6, 0.40105, 0.1, 0.23975])

n = len(x) - 1
h = np.array([x[k + 1] - x[k] for k in range(n)])
system = np.identity(len(x))
vector = np.zeros(len(x))

for i in range(1, n):
    system[i][i - 1] = h[i - 1]
    system[i][i] = 2 * (h[i - 1] + h[i])
    system[i][i + 1] = h[i]
    vector[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

c = np.array(np.linalg.solve(system, vector))
b = np.zeros(len(x))
d = np.zeros(len(x))
for i in range(1, len(x)):
    b[i] = (y[i] - y[i - 1]) / h[i - 1] + (2 * c[i] + c[i - 1]) * h[i - 1] / 3
    d[i] = (c[i] - c[i - 1]) / h[i - 1] / 3

# S[i](x) = a[i] + b[i](x - x[i]) + c[i](x - x[i])**2 + d[i](x - x[i])**3
def S(value):
    if value == x[n]:
        return y[n]
    i = binSearch(x, value)
    return y[i] \
           + b[i]*(value - x[i]) \
           + c[i]*(value - x[i])**2 \
           + d[i]*(value - x[i])**3

xRange = np.linspace(x[0], x[n], 100)
plt.plot(x, y, 'ro')
plt.plot(xRange, [S(value) for value in xRange])
plt.show()