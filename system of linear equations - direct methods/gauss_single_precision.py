import numpy as np
import sys


def matrixConditioning(matrix):
    inv_matrix = np.linalg.inv(matrix)
    return np.linalg.norm(matrix, np.inf) * np.linalg.norm(inv_matrix, np.inf)


#   x -> n - elem vector with elements = 1
#   a -> n x n matrix with a1j = 1 and aij = 1 / (i + j - 1)
def buildMatricesSinglePrecision1(size):
    a = np.zeros((size, size))
    x = np.zeros(size)

    for i in range(size):
        x[i] = 1
        for j in range(size):
            if i == 0:
                #   changes in matrix elements formula result from python array indexing (0 to size - 1, not 1 to size)
                a[i][j] = 1
            else:
                a[i][j] = np.single(1 / (i + j + 1))
    b = np.single(np.matmul(a, x))
    return a, x, b


#   x -> n - elem vector with elements = 1
#   a -> n x n matrix with aij = 2i / j for j >= 1 and aij = aji for j < i (i, j = 1,...,n)
def buildMatricesSinglePrecision2(size):
    a = np.zeros((size, size))
    x = np.zeros(size)

    for i in range(size):
        x[i] = 1
        for j in range(size):
            if j >= i:
                a[i][j] = np.single(2 * (i + 1) / (j + 1))
            else:
                a[i][j] = np.single(a[j][i])
    b = np.single(np.matmul(a, x))
    return a, x, b


def forwardEliminationSinglePrecision(augmented_a, size):
    for i in range(size):
        if augmented_a[i][i] == 0:
            k = i + 1
            while k < size:
                if augmented_a[k][i] != 0:
                    augmented_a[[i, k]] = augmented_a[[k, i]]
                    break

                k += 1

            if augmented_a[i][i] == 0:
                sys.exit('Division by 0')

        for j in range(i + 1, size):
            ratio = np.single(augmented_a[j][i] / augmented_a[i][i])

            for k in range(size + 1):
                augmented_a[j][k] = np.single(augmented_a[j][k] - np.single(ratio * augmented_a[i][k]))

    if augmented_a[-1][-1] == 0:
        sys.exit('No unique solution exists')


def backwardSubstitutionSinglePrecision(augmented_a, size):
    x = np.zeros(size)
    x[-1] = np.single(augmented_a[-1][-1] / augmented_a[-1][size - 1])

    for i in range(size - 2, -1, -1):
        x[i] = np.single(augmented_a[i][-1])

        for j in range(i + 1, size):
            x[i] = np.single(x[i] - np.single(augmented_a[i][j] * x[j]))

        x[i] = np.single(x[i] / augmented_a[i][i])

    return x


def gaussianElimination(size, matrix):
    if matrix == 1:
        a, x, b = buildMatricesSinglePrecision1(size)
    elif matrix == 2:
        a, x, b = buildMatricesSinglePrecision2(size)

    cond = matrixConditioning(a)

    report.write("\tMatrix conditioning: %10.3E\n" % cond)
    print("Matrix conditioning: %10.3E" % cond)

    augmented_a = np.insert(a, size, b, axis=1)

    forwardEliminationSinglePrecision(augmented_a, size)

    calculated_x = backwardSubstitutionSinglePrecision(augmented_a, size)

    report.write("\tEuclidean norm: %10.3E\n" % np.linalg.norm(calculated_x - x))
    print("Euclidean norm: %10.3E" % np.linalg.norm(calculated_x - x))
    report.write("\tMax norm: %10.3E\n" % np.linalg.norm(calculated_x - x, np.inf))
    print("Max norm: %10.3E" % np.linalg.norm(calculated_x - x, np.inf))


if __name__ == '__main__':
    n = int(sys.argv[1])
    matrix = int(sys.argv[2])
    if matrix not in [1, 2]:
        sys.exit("Wrong matrix type")
    if n <= 0:
        sys.exit("Wrong matrix size")
    report = open("gauss_single_precision.txt", "a")
    report.write("\nGauss single precision for matrix size n = %d and mode = %d\n" % (n, matrix))
    gaussianElimination(n, matrix)
