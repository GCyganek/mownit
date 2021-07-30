import numpy as np
import sys
import time

def matrixConditioning(matrix):
    inv_matrix = np.linalg.inv(matrix)
    return np.linalg.norm(matrix, np.inf) * np.linalg.norm(inv_matrix, np.inf)


#   aii = -5i - 5
#   ai(i+1) = i
#   ai(i-1) = 5 / i for i > 1
#   aij = 0 for j < i - 1 and j > i + 1
def buildMatricesDoublePrecision(size):
    a = np.zeros((size, size))
    x = np.zeros(size)
    bottom = np.zeros(size - 1)
    middle = np.zeros(size)
    top = np.zeros(size - 1)

    for i in range(size):
        x[i] = 1
        a[i][i] = middle[i] = -5 * (i + 1) - 5
        if i - 1 >= 0:
            a[i][i - 1] = bottom[i - 1] = 5 / (i + 1)
        if i + 1 < size:
            a[i][i + 1] = top[i] = (i + 1)
    b = np.matmul(a, x)
    return a, bottom, middle, top, x, b


def thomasAlgorithmWithTime(size):
    a, bottom, middle, top, x, b = buildMatricesDoublePrecision(size)

    cond = matrixConditioning(a)

    report.write("\tMatrix conditioning: %10.3E\n" % cond)
    print("Matrix conditioning: %10.3E" % cond)

    start = time.perf_counter()

    # thomas algorithm starts here...
    for i in range(1, size):
        multiplier = bottom[i - 1] / middle[i - 1]  # its bottom in i row, but we have array without bottom element in first row which doesnt exist
        middle[i] -= (multiplier * top[i - 1])
        b[i] -= (multiplier * b[i - 1])

    calculated_x = middle
    calculated_x[-1] = b[-1] / middle[-1]

    for i in range(size - 2, -1, -1):
        calculated_x[i] = (b[i] - top[i] * calculated_x[i + 1]) / middle[i]
    # ... and ends here

    end = time.perf_counter()

    report.write("\tEuclidean norm: %10.3E\n" % np.linalg.norm(calculated_x - x))
    print("Euclidean norm: %10.3E" % np.linalg.norm(calculated_x - x))
    report.write("\tMax norm: %10.3E\n" % np.linalg.norm(calculated_x - x, np.inf))
    print("Max norm: %10.3E" % np.linalg.norm(calculated_x - x, np.inf))

    report.write("\tTime used to receive the solution: " + str(end - start))
    print("Time used to receive the solution: " + str(end - start))


if __name__ == '__main__':
    n = int(sys.argv[1])
    if n <= 0:
        sys.exit("Wrong matrix size")
    report = open("thomas_algorithm.txt", "a")
    report.write("\nGauss double precision for matrix size n = %d\n" % n)
    thomasAlgorithmWithTime(n)

