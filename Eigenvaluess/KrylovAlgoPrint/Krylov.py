import numpy as np

def myPrint(obj, title = ""):
    print(title)
    try:
        obj.shape[1]
        for row in obj:
            for elem in row:
                print("%10.4f" % elem, end=' ')
            print()
    except IndexError:
        for elem in obj:
            print("%10.4f" % elem, end = ' ')
        print()
    print()

with open("input.txt", "r") as file:
    n = int(file.readline())
    matrix = ""
    for i in range(n):
        matrix = matrix + file.readline().strip() + " "
    matrix = np.array([float(word) for word in matrix.split()]).reshape(n, n)

c = np.zeros((n + 1, n), float)
c[0] = np.identity(n, float)[0, ]
for i in range(1, n + 1):
    c[i] = np.dot(matrix, c[i - 1])

matrixSystem = np.transpose(np.concatenate(c[n-1::-1]).reshape(n, n))
solution = np.linalg.solve(matrixSystem, c[n])
polynomial = np.concatenate((np.array([1]), solution * -1))
eigenvalues = np.roots(polynomial)

myPrint(c.transpose(), "C-vectors: ")
myPrint(matrixSystem, "Matrix-system: ")
myPrint(c[n], "Vector of the system: ")
myPrint(polynomial, "Coefficients of polynomial: ")
myPrint(eigenvalues, "Eigenvalues of matrix: ")

for eigenvalue in eigenvalues:
    beta = np.zeros(n, float)
    beta[0] = 1
    eigenvector = np.zeros(n, float)
    for i in range(1, n):
        beta[i] = beta[i - 1] * eigenvalue + polynomial[i]
    myPrint(beta, "Beta for %-.4f" % eigenvalue)
    for i in range(0, n):
        eigenvector += beta[i] * matrixSystem[0:n, i]
    myPrint(eigenvector, "Eigenvector for %-.4f" % eigenvalue)
