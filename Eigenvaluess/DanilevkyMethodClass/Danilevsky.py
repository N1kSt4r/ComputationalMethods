import numpy as np


class DanilevskyMethod:
    def __init__(self, matrix):
        self.n = len(matrix)
        self.matrix = np.copy(matrix)
        self.recalculate()

    def recalculate(self):
        self.__transform_matrix()
        self.__calculate_values()
        self.__calculate_vectors()

    def __transform_matrix(self):
        self.resultingS = np.identity(self.n)
        for j in range(self.n - 1, 0, -1):
            S = np.identity(self.n)
            S[j - 1] = -self.matrix[j] / self.matrix[j, j - 1]
            S[j - 1, j - 1] /= -self.matrix[j, j - 1]
            self.resultingS = np.dot(self.resultingS, S)
            self.matrix = np.dot(np.linalg.inv(S), self.matrix)
            self.matrix = np.dot(self.matrix, S)

    def __calculate_values(self):
        self.eigenvalues = np.roots(np.concatenate((np.ones(1), self.matrix[0] * -1)))

    def __calculate_vectors(self):
        self.eigenvectors = np.zeros((self.n, self.n), float)
        for i, eigenvalue in enumerate(self.eigenvalues):
            eigenvector = np.ones(self.n)
            for j in range(self.n - 1, 0, -1):
                eigenvector[j - 1] = eigenvalue * eigenvector[j]
            eigenvector = np.dot(self.resultingS, eigenvector)
            self.eigenvectors[i] = eigenvector
