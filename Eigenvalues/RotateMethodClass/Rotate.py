import numpy as np


class RotateMethod:
    def __init__(self, matrix, eps = 0.1 ** 7):
        self.matrix = matrix
        self.eps = eps
        self.n = len(matrix)
        self.eigenvectors = np.identity(self.n)
        self.calculate()

    def __zerodiagonal(self):
        temp = np.copy(self.matrix)
        for i in range(self.n):
            temp[i, i] = 0
        return temp

    def __Givens(self, cos, sin, i, j):
        Givens = np.identity(self.n)
        Givens[i, i] = Givens[j, j] = cos
        Givens[i, j] = Givens[j, i] = sin
        Givens[i, j] *= -1
        return Givens

    def calculate(self):
        self.eigenvectors = np.identity(self.n)
        while np.linalg.norm(self.__zerodiagonal()) > self.eps:
            posmax = np.argmax(np.abs(self.__zerodiagonal()))
            i, j = posmax // self.n, posmax % self.n
            temp = 2 * self.matrix[i, j] / (self.matrix[i, i] - self.matrix[j, j])
            temp = (1 + temp ** 2) ** -0.5
            cos = ((1 + temp) / 2) ** 0.5
            sin = ((1 - temp) / 2) ** 0.5
            givens = self.__Givens(cos, sin, i ,j)
            self.eigenvectors = np.dot(self.eigenvectors, givens)
            self.matrix = np.dot(givens.T, np.dot(self.matrix, givens))
        self.eigenvalues = np.diag(self.matrix)
        self.eigenvectors = self.eigenvectors.T
