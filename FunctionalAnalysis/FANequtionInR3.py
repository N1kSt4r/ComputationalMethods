import scipy as sp
from scipy.misc import derivative
import numpy as np
import sympy
import matplotlib.pyplot as plt

# Тема 3, Задание 1, Вариант 21 (7)
# 1.7
#  0.9x1 + 0.1x2 + 0.1x3 = 1
#  0.2x1 + 1.1x2 + 0.2x3 = 1
#  0.1x1 + 0.2x2 + 1.2x3 = 0

xSymbol = sympy.Symbol('x')
eps = 1e-3
A = np.array([
    [0.9, 0.1, 0.1],
    [0.2, 1.1, 0.2],
    [0.1, 0.2, 1.2]
])
b = np.array([1, 1, 0])

B = -A + np.identity(len(A))
y = -b
x = np.zeros(len(A))

print("A:\n", A)
print("b: ", b)
print("B:\n", B)
print("y: ", y)
normB = np.linalg.norm(B)
print("||B|| =", normB)
nextApproximation = lambda B, x, y: np.dot(B, x) - y
estimate = int(np.log(eps * (1 - normB) / np.linalg.norm(x - nextApproximation(B, x, y))) / np.log(normB)) + 1
print("Примерная оценка:", estimate)

for i in range(estimate):
    print((" %+.4f" * len(x)) % tuple(value for value in x))
    x = nextApproximation(B, x, y)

print("Решение:", x)
print("F(x) =", np.dot(A, x) - b)
