import numpy as np


class KrylovMethod:
    def __init__(self, matrix):
        self.__n = len(matrix)
        self.__matrix = np.copy(matrix)
        self.__values_called = False
        self.__vector_called = False

    def __calculate_c(self):
        c = np.zeros((self.__n + 1, self.__n), float)
        c[0] = np.identity(n, float)[0,]
        for i in range(1, self.__n + 1):
            c[i] = np.dot(self.__matrix, c[i - 1])
        return c

    def __calculate_beta(self, eigenvalue):
        beta = np.zeros(self.__n, float)
        beta[0] = 1
        for i in range(1, self.__n):
            beta[i] = beta[i - 1] * eigenvalue + self.__polynomial[i]
        return beta

    def get_eigenvalues(self):
        if self.__values_called:
            return np.copy(self.__eigenvalues)
        
        c = self.__calculate_c()
        self.__matrixSystem = np.transpose(np.concatenate(c[self.__n - 1::-1]).reshape(self.__n, self.__n))
        solution = np.linalg.solve(self.__matrixSystem, c[self.__n])
        self.__polynomial = np.concatenate((np.array([1]), solution * -1))
        self.__eigenvalues = np.roots(self.__polynomial)
        
        self.__values_called = True
        return np.copy(self.__eigenvalues)

    def get_eigenvectors(self):
        if self.__vector_called:
            return np.copy(self.__eigenvectors)
        if not self.__values_called:
            self.get_eigenvalues()
        self.__eigenvectors = np.zeros(self.__n, float)
        
        for eigenvalue in self.__eigenvalues:
            beta = self.__calculate_beta(eigenvalue)
            eigenvector = np.zeros(self.__n, float)
            for i in range(0, self.__n):
                eigenvector += beta[i] * self.__matrixSystem[0:self.__n, i]
            self.__eigenvectors = np.concatenate((self.__eigenvectors, eigenvector))
            
        self.__vector_called = True
        return np.copy(np.reshape(self.__eigenvectors[self.__n:], (self.__n, self.__n)))
