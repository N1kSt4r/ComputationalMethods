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
    myPrint(matrix, "Original matrix: ")

eps = 0.1 ** 7

def zerodiagonal(matrix):
    temp = np.copy(matrix)
    for i in range(len(temp)):
        temp[i, i] = 0
    return temp

GivensRes = np.identity(n)

while np.linalg.norm(zerodiagonal(matrix)) > eps:
    posmax = np.argmax(np.abs(zerodiagonal(matrix)))
    i, j = posmax // n, posmax % n
    temp = 2 * matrix[i, j] / (matrix[i, i] - matrix[j, j])
    temp = (1 + temp ** 2) ** -0.5
    cos = ((1 + temp) / 2) ** 0.5
    sin = ((1 - temp) / 2) ** 0.5
    Givens = np.identity(n)
    Givens[i, i] = Givens[j, j] = cos
    Givens[i, j] = Givens[j, i] = sin
    Givens[i, j] *= -1
    GivensRes = np.dot(GivensRes, Givens)
    matrix = np.dot(matrix, Givens)
    matrix = np.dot(Givens.T, matrix)

GivensRes = GivensRes.T
eigenvalue = np.diag(matrix)

myPrint(matrix, "Transformed matrix: ")
myPrint(eigenvalue, "Eigenvalues: ")
for i in range(n):
    myPrint(GivensRes[i], "Eigenvector for %.4f: " % eigenvalue[i])
