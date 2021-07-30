import math
import numpy
import matplotlib.pyplot as plt

def functionValue(x):
    return (x - 1) * math.exp(-10 * x) + math.pow(x, 15)

def functionValues(nodes):
    result = []
    for xk in nodes:
        yk = functionValue(xk)
        result.append(yk)
    return result

def drawPlot():
    nodes = numpy.linspace(-1, 0.7, 1000)
    f_values = functionValues(nodes)

    plt.axhline(y=0, color='k', linewidth=1)
    plt.axvline(x=0, color='k', linewidth=1)
    plt.plot(nodes, f_values, 'm-', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Zadana funkcja na przedziale [-1, 0.7]")
    # plt.ylim(-0.5, 0.5)
    # plt.xlim(-0.1, 0.8)
    plt.show()


drawPlot()
