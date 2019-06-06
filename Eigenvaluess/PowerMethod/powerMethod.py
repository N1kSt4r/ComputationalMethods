import numpy as np

def myPrint(obj, title = ""):
    print(title)
    try:
        obj.shape[1]
        for row in obj:
            for elem in row:
                print("%10.4f" % elem, end = ' ')
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

eigenvectorTemp = np.dot(matrix, np.ones(n))
eigenvector = np.dot(matrix, eigenvectorTemp)

ratio = eigenvector / eigenvectorTemp

matrixPow2 = np.dot(matrix, matrix)
while True:
    eigenvector = eigenvector / max(abs(eigenvector))
    eigenvectorTemp = eigenvector
    eigenvector = np.dot(matrixPow2, eigenvector)
    ratioTemp = ratio
    ratio = eigenvector / eigenvectorTemp
    if np.linalg.norm(ratio - ratioTemp) < eps:
        break

print("Maximal eigenvalue: \n%10.4f\n" % ratio[0] ** 0.5)
myPrint(eigenvector / max(abs(eigenvector)), "Eigenvector: ")
#eigenvector correct only if matrix hasn't negative of same maximal eigenvalue

