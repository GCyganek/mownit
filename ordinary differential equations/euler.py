import matplotlib.pyplot as plt
import sys
from utils import *

def eulerMethod(x, k, y0):
    h = x[1] - x[0]

    y = np.zeros(k)
    y[0] = y0
    real_y = np.zeros(k)
    real_y[0] = realFunction(x[0])

    max_err = 0

    for i in range(1, k):
        y[i] = y[i - 1] + function(x[i - 1], y[i - 1]) * h
        real_y[i] = realFunction(x[i])
        max_err = max(abs(y[i] - real_y[i]), max_err)

    print("Euclidean norm: %10.3E" % sqrError(y, real_y))
    print("Max norm: %10.3E" % max_err)

    return y


if __name__ == '__main__':
    k = int(sys.argv[1])
    if k < 2:
        print("Higher k value needed (k >= 2)")
        sys.exit()
    x = np.linspace(x0, xk, k)
    y0 = realFunction(x0)
    y = eulerMethod(x, k, y0)

    plt_x = np.linspace(x0, xk, 1000)
    plt_y = realFunctionValues(plt_x, 1000)

    plt.plot(plt_x, plt_y, 'm--', label='Real function')
    plt.plot(x, y, 'b-', label='Euler method result')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Euler method for %d points" % k)
    plt.legend(loc='lower center', fontsize='x-small')
    plt.savefig("euler/Euler%d" % k)
