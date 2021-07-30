import numpy as np
import sys
import time


def matrixConditioning(matrix):
    inv_matrix = np.linalg.inv(matrix)
    return np.linalg.norm(matrix, np.inf) * np.linalg.norm(inv_matrix, np.inf)


#   x -> n - elem vector with elements = 1
#   a -> n x n matrix with a1j = 1 and aij = 1 / (i + j - 1)
def buildMatricesDoublePrecision1(size):
    a = np.zeros((size, size))
    x = np.zeros(size)

    for i in range(size):
        x[i] = 1
        for j in range(size):
            if i == 0:
                #   changes in matrix elements formula result from python array indexing (0 to size - 1, not 1 to size)
                a[i][j] = 1
            else:
                a[i][j] = 1 / (i + j + 1)
    b = np.matmul(a, x)
    return a, x, b


#   x -> n - elem vector with elements = 1
#   a -> n x n matrix with aij = 2i / j for j >= 1 and aij = aji for j < i (i, j = 1,...,n)
def buildMatricesDoublePrecision2(size):
    a = np.zeros((size, size))
    x = np.zeros(size)

    for i in range(size):
        x[i] = 1
        for j in range(size):
            if j >= i:
                a[i][j] = 2 * (i + 1) / (j + 1)
            else:
                a[i][j] = a[j][i]
    b = np.matmul(a, x)
    return a, x, b


#   aii = -5i - 5
#   ai(i+1) = i
#   ai(i-1) = 5 / i for i > 1
#   aij = 0 for j < i - 1 and j > i + 1
def buildMatricesDoublePrecision3(size):
    a = np.zeros((size, size))
    x = np.zeros(size)

    for i in range(size):
        x[i] = 1
        a[i][i] = -5 * (i + 1) - 5
        if i - 1 >= 0:
            a[i][i - 1] = 5 / (i + 1)
        if i + 1 < size:
            a[i][i + 1] = (i + 1)
    b = np.matmul(a, x)
    return a, x, b


def forwardEliminationDoublePrecision(augmented_a, size):
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
            ratio = augmented_a[j][i] / augmented_a[i][i]

            for k in range(size + 1):
                augmented_a[j][k] = augmented_a[j][k] - ratio * augmented_a[i][k]

    if augmented_a[-1][-1] == 0:
        sys.exit('No unique solution exists')


def backwardSubstitutionDoublePrecision(augmented_a, size):
    x = np.zeros(size)
    x[-1] = augmented_a[-1][-1] / augmented_a[-1][size - 1]

    for i in range(size - 2, -1, -1):
        x[i] = augmented_a[i][-1]

        for j in range(i + 1, size):
            x[i] -= augmented_a[i][j] * x[j]

        x[i] /= augmented_a[i][i]

    return x


def gaussianElimination(size, matrix):
    if matrix == 1:
        a, x, b = buildMatricesDoublePrecision1(size)
    elif matrix == 2:
        a, x, b = buildMatricesDoublePrecision2(size)

    cond = matrixConditioning(a)

    report.write("\tMatrix conditioning: %10.3E\n" % cond)
    print("Matrix conditioning: %10.3E" % cond)

    augmented_a = np.insert(a, size, b, axis=1)

    forwardEliminationDoublePrecision(augmented_a, size)

    calculated_x = backwardSubstitutionDoublePrecision(augmented_a, size)

    report.write("\tEuclidean norm: %10.3E\n" % np.linalg.norm(calculated_x - x))
    print("Euclidean norm: %10.3E" % np.linalg.norm(calculated_x - x))
    report.write("\tMax norm: %10.3E\n" % np.linalg.norm(calculated_x - x, np.inf))
    print("Max norm: %10.3E" % np.linalg.norm(calculated_x - x, np.inf))


def gaussianEliminationWithTime(size):
    a, x, b = buildMatricesDoublePrecision3(size)

    cond = matrixConditioning(a)

    report.write("\tMatrix conditioning: %10.3E\n" % cond)
    print("Matrix conditioning: %10.3E" % cond)

    start = time.perf_counter()

    augmented_a = np.insert(a, size, b, axis=1)

    forwardEliminationDoublePrecision(augmented_a, size)

    calculated_x = backwardSubstitutionDoublePrecision(augmented_a, size)

    end = time.perf_counter()

    report.write("\tEuclidean norm: %10.3E\n" % np.linalg.norm(calculated_x - x))
    print("Euclidean norm: %10.3E" % np.linalg.norm(calculated_x - x))
    report.write("\tMax norm: %10.3E\n" % np.linalg.norm(calculated_x - x, np.inf))
    print("Max norm: %10.3E" % np.linalg.norm(calculated_x - x, np.inf))

    report.write("\tTime used to receive the solution: " + str(end - start))
    print("Time used to receive the solution: " + str(end - start))


if __name__ == '__main__':
    n = int(sys.argv[1])
    matrix = int(sys.argv[2])
    if matrix not in [1, 2, 3]:
        sys.exit("Wrong matrix type")
    if n <= 0:
        sys.exit("Wrong matrix size")
    report = open("gauss_double_precision.txt", "a")
    report.write("\nGauss double precision for matrix size n = %d and mode = %d\n" % (n, matrix))
    if matrix == 3:
        gaussianEliminationWithTime(n)
    else:
        gaussianElimination(n, matrix)
