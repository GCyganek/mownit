import matplotlib.pyplot as plt
import sys
from utils import *

def rungeKuttaMethod(x, k, y0, degree):
    h = x[1] - x[0]

    y = np.zeros(k)
    y[0] = y0
    real_y = np.zeros(k)
    real_y[0] = realFunction(x[0])

    max_err = 0

    def rungeKuttaFunction(index):
        k1 = h * function(x[index], y[index])
        if degree == 1:
            return k1

        k2 = h * function(x[index] + h / 2, y[index] + k1 / 2)
        if degree == 2:
            return k2

        k3 = h * function(x[index] + h / 2, y[index] + k2 / 2)
        if degree == 3:
            return k3

        k4 = h * function(x[index] + h, y[index] + k3)
        if degree == 4:
            return (k1 + 2 * k2 + 2 * k3 + k4) / 6

    for i in range(1, k):
        y[i] = y[i - 1] + rungeKuttaFunction(i - 1)
        real_y[i] = realFunction(x[i])
        max_err = max(abs(y[i] - real_y[i]), max_err)

    print("Euclidean norm: %10.3E" % sqrError(y, real_y))
    print("Max norm: %10.3E" % max_err)

    return y


if __name__ == '__main__':
    k = int(sys.argv[1])
    degree = int(sys.argv[2])
    if k < 2:
        print("Higher k value needed (k >= 2)")
        sys.exit()
    x = np.linspace(x0, xk, k)
    y0 = realFunction(x0)
    y = rungeKuttaMethod(x, k, y0, degree)

    plt_x = np.linspace(x0, xk, 1000)
    plt_y = realFunctionValues(plt_x, 1000)

    plt.plot(plt_x, plt_y, 'm--', label='Real function')
    plt.plot(x, y, 'b-', label='Runge-Kutta method result')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Runge-Kutta method for %d points, %d degree" % (k, degree))
    plt.legend(loc='lower center', fontsize='x-small')
    # plt.show()
    plt.savefig("runge-kutta/RungeKutta%d-deg%d" % (k, degree))
