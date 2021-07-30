import sys
from lab_utils import *
import time

# spectral radius of D^(-1) * (L + U)
def spectralRadiusJacobi(a):
    size = len(a)
    d_inv = np.zeros((size, size))

    for i in range(size):
        d_inv[i][i] = 1 / a[i][i]
        a[i][i] = 0

    report.write("\tspectral radius of iteration matrix: %10.3E\n" % max(abs(np.linalg.eigvals(np.matmul(d_inv, a)))))
    print("\tspectral radius of iteration matrix: %10.3E" % max(abs(np.linalg.eigvals(np.matmul(d_inv, a)))))


def jacobiMethod(a, x0, b, ro, condition, x):
    start = time.perf_counter()
    step = 0
    x1 = np.zeros(len(x0))
    while True:
        for i in range(len(x0)):
            multiplier = 1 / a[i][i]
            sum = 0
            for j in range(len(x0)):
                if i != j:
                    sum -= a[i][j] * x0[j]
            sum *= multiplier
            sum += multiplier * b[i]
            x1[i] = sum

        step += 1
        if condition == 1:
            if condition1(x0, x1, ro):
                break
        elif condition == 2:
            if condition2(a, x1, b, ro):
                break

        x0 = x1.copy()
    end = time.perf_counter()
    err = np.linalg.norm(x1 - x)
    report.write("\titerations = %d, error = %10.3E\n" % (step, err))
    print("\titerations = %d, error = %10.3E" % (step, err))
    report.write("\tTime used to receive the solution: " + str(end - start))
    print("\ttime used to receive the solution: " + str(end - start))


if __name__ == '__main__':
    report = open("jacobi_method.txt", "a")
    n = int(sys.argv[1])
    cond = int(sys.argv[2])
    ro = float(sys.argv[3])
    a, x, b = buildMatrices(n)

    if cond == 1:
        cond_str = "||x^(i+1) - x^(i)|| < ro"
    elif cond == 2:
        cond_str = "||Ax^(i) - b|| < ro"
    else:
        sys.exit()

    if n < 1:
        sys.exit()

    report.write("\nJacobi method for n = %d, ro = %10.3E, condition: %s and x0 with 0 values\n"
                 % (n, ro, cond_str))
    print("Jacobi method for n = %d, ro = %10.3E, condition: %s and x0 with 0 values"
          % (n, ro, cond_str))

    x0 = np.zeros(n)

    if len(sys.argv) == 5 and sys.argv[4] == "print":
        print("\tInitial vector x0: %s" % x0)
        report.write("\tInitial vector x0: %s\n" % x0)

    jacobiMethod(a, x0, b, ro, cond, x)
    spectralRadiusJacobi(a)
