import sys
from lab_utils import *
import time

# spectral radius of (D - omega * L)^(-1) * [(1 - omega) * D + omega * U]
def spectralRadiusSOR(a):
    size = len(a)
    d = np.zeros((size, size))
    u = np.zeros((size, size))
    l = np.zeros((size, size))

    for i in range(size):
        d[i][i] = a[i][i]
        for j in range(size):
            if i < j:
                u[i][j] = a[i][j]
            if i > j:
                l[i][j] = a[i][j]

    inv = np.linalg.inv(d - omega * l)
    result = np.matmul(inv, ((1 - omega) * d + omega * u))

    report.write("\tspectral radius of iteration matrix: %10.3E\n" % max(abs(np.linalg.eigvals(result))))
    print("\tspectral radius of iteration matrix: %10.3E" % max(abs(np.linalg.eigvals(result))))


def sorMethod(a, x0, b, ro, condition, x, omega):
    start = time.perf_counter()
    step = 0
    while True:
        x_to_compare = x0.copy()
        for i in range(len(x0)):
            sum = 0
            for j in range(len(x0)):
                if i != j:
                    sum += a[i][j] * x0[j]
            x0[i] = (1 - omega) * x0[i] + (omega / a[i][i]) * (b[i] - sum)

        step += 1
        if condition == 1:
            if condition1(x_to_compare, x0, ro):
                break
        elif condition == 2:
            if condition2(a, x0, b, ro):
                break

    end = time.perf_counter()
    err = np.linalg.norm(x0 - x)
    report.write("\titerations = %d, error = %10.3E\n" % (step, err))
    print("\titerations = %d, error = %10.3E" % (step, err))
    report.write("\tTime used to receive the solution: " + str(end - start))
    print("\ttime used to receive the solution: " + str(end - start))


if __name__ == '__main__':
    report = open("sor_method.txt", "a")
    n = int(sys.argv[1])
    cond = int(sys.argv[2])
    ro = float(sys.argv[3])
    omega = float(sys.argv[4])
    a, x, b = buildMatrices(n)

    if cond == 1:
        cond_str = "||x^(i+1) - x^(i)|| < ro"
    elif cond == 2:
        cond_str = "||Ax^(i) - b|| < ro"
    else:
        sys.exit()

    if n < 1 or omega <= 0 or omega >= 2:
        sys.exit()

    report.write("\nSor method for n = %d, ro = %10.3E, condition: %s, omega = %d and x0 with 0 values\n "
                 % (n, ro, cond_str, omega))
    print("Sor method for n = %d, ro = %10.3E, condition: %s, omega = %d and x0 with 0 values"
          % (n, ro, cond_str, omega))

    x0 = np.zeros(n)

    if len(sys.argv) == 6 and sys.argv[5] == "print":
        print("\tInitial vector x0: %s" % x0)
        report.write("\tInitial vector x0: %s\n" % x0)

    sorMethod(a, x0, b, ro, cond, x, omega)
    spectralRadiusSOR(a)
