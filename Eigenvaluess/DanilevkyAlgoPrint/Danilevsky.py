import numpy as np
import mpmath

mpmath.mp.dps = 30

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
    matrix = np.array([mpmath.mpf(word) for word in matrix.split()]).reshape(n, n)
    myPrint(matrix, "Original matrix: ")

regularCase = True
resultingS = np.identity(n)
indexCeil = [0]
for j in range(n - 1, 0, -1):
    S = np.identity(n)
    i = j - 1
    while i >= 0 and abs(matrix[j, i]) < 1e-30:
        i = i - 1
    if i is not -1:
        if i is not j - 1:
            tempS = np.identity(n)
            tempS[i], tempS[j - 1] = np.copy(tempS[j - 1]), np.copy(tempS[i])
            myPrint(tempS, "Additional transform matrix (%d): " % j)
            matrix = np.dot(tempS.T, matrix)
            matrix = np.dot(matrix, tempS)
            resultingS = np.dot(resultingS, tempS)
        S[j - 1] = -matrix[j] / matrix[j, j - 1]
        S[j - 1, j - 1] /= -matrix[j, j - 1]
        myPrint(S, "Transform matrix (%d): " % j)
        resultingS = np.dot(resultingS, S)
        matrix = np.dot(np.linalg.inv(S), matrix)
        matrix = np.dot(matrix, S)
        myPrint(matrix, "Matrix: ")
    else:
        indexCeil.insert(len(indexCeil), j)
        regularCase = False


if regularCase:
    polynomial = np.concatenate((np.array([1]), matrix[0] * -1))
    eigenvalues = np.roots(polynomial)

    myPrint(resultingS, "Resulting transform matrix: ")
    myPrint(matrix, "Matrix after transformations: ")
    myPrint(polynomial, "Coefficients of polynomial: ")
    myPrint(eigenvalues, "Eigenvalues: ")
    for eigenvalue in eigenvalues:
        eigenvector = np.ones(n)
        for j in range(n - 1, 0, -1):
            eigenvector[j - 1] = eigenvalue * eigenvector[j]
        eigenvector = np.dot(resultingS, eigenvector)
        myPrint(eigenvector, "Eigenvector for %.4f: " % eigenvalue)
else:
    print("Irregular case\n")
    indexCeil.insert(len(indexCeil), n)
    for i in range(len(indexCeil) - 1):
        temp = [-1]
        for j in range(indexCeil[i], indexCeil[i + 1]):
            temp.insert(len(temp), matrix[indexCeil[i], j])
        myPrint(np.roots(np.array(temp)), "Values of %d ceil:" % (i + 1))
