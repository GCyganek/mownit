import sys
import math
import numpy
import matplotlib.pyplot as plot


# Lagrange and Newton Interpolation for f(x) = exp(-2sin(2x)) + 2cos(2x), [-pi, 2pi]


# xk = 1/2(a + b) + 1/2(b - a)cos((2k - 1)pi/2n), k=1...n for [a, b]
# for [-pi, 2pi] it will be
# xk = 1/2(-pi + 2pi) + 1/2(2p + pi)cos((2k - 1)pi/2n), k=1...n for [-pi, 2pi]
# xk = pi/2 + (3pi/2)cos((2k-1)pi/2n), k=1...n for [-pi, 2pi]
def chebyshevNodesFromMinusPiTo2Pi(nodes_number):
    result = []
    for k in range(nodes_number, 0, -1):
        xk = math.pi / 2 + (3 * math.pi / 2) * math.cos(((2 * k - 1) * math.pi) / (2 * nodes_number))
        result.append(xk)
    return result


def equallySpacedNodesFromMinusPiTo2Pi(nodes_number):
    result = []
    space_between = (3 * math.pi) / (nodes_number - 1)
    for k in range(0, nodes_number):
        xk = -math.pi + space_between * k
        result.append(xk)
    return result


def valueOfFunctionInPoint(x):
    return math.exp(-2 * math.sin(2 * x)) + 2 * math.cos(2 * x)


def valuesOfFunctionInNodes(nodes):
    result = []
    for xk in nodes:
        yk = valueOfFunctionInPoint(xk)
        result.append(yk)
    return result


def valueOfDerivativeInPoint(x):
    return -4 * math.exp(-2 * math.sin(2 * x)) * math.cos(2 * x) - 4 * math.sin(2 * x)


def valuesOfDerivativeInNodes(nodes):
    result = []
    for xk in nodes:
        yk = valueOfDerivativeInPoint(xk)
        result.append(yk)
    return result


def lagrangeInterpolation(nodes):
    x = numpy.array(nodes, float)
    y = numpy.array(valuesOfFunctionInNodes(x), float)
    xplt = numpy.linspace(-math.pi, 2 * math.pi, num=900)
    yplt = numpy.array([], float)
    max_err = 0

    for xp in xplt:
        yp = 0
        y_fun = valueOfFunctionInPoint(xp)

        for xi, yi in zip(x, y):
            yp += yi * numpy.product((xp - x[x != xi]) / (xi - x[x != xi]))
        yplt = numpy.append(yplt, yp)
        max_err = max(abs(y_fun - yp), max_err)

    y_fun = valuesOfFunctionInNodes(xplt)
    sqr_err = numpy.square(numpy.subtract(y_fun, yplt)).mean()
    sqr_err = math.sqrt(sqr_err)
    return x, y, xplt, yplt, max_err, sqr_err

def hermiteLagrangeInterpolation(nodes):
    x = numpy.array(nodes, float)
    y = numpy.array(valuesOfFunctionInNodes(x), float)
    xplt = numpy.linspace(-math.pi, 2 * math.pi, num=900)
    yplt = numpy.array([], float)

    for xp in xplt:
        lagrange = 0

        for xi in x:
            lagrange += numpy.product((xp - x[x != xi]) / (xi - x[x != xi]))

        yp = 0

        for xi, yi in zip(x, y):
            yp += yi * ((1 - 2 * (xp - xi) * lagrange * numpy.sum(1 / (x[x != xi] - xi))) * (lagrange ** 2))
            yp += valueOfDerivativeInPoint(xi) * numpy.power(lagrange, 2) * (xp - xi)
        yplt = numpy.append(yplt, yp)

    return x, y, xplt, yplt


def lagrangeInterpolationChebyshev(nodes_number):
    nodes = chebyshevNodesFromMinusPiTo2Pi(nodes_number)
    (x, y, xplt, yplt, max_err, sqr_err) = lagrangeInterpolation(nodes)

    plot.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Lagrange`a używająca ' + str(nodes_number) + ' węzłów Czebyszewa')
    plot.legend(('Węzły', 'wielomian Lagrange`a', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/lagrange_chebyshev_plot_%d' % nodes_number)
    plot.close()

    report.write("Lagrange %d wezlow Czebyszewa:" % nodes_number + "         max_error =%10.3E" % max_err +
                 " sqr_err =%10.3E\n" % sqr_err)


def lagrangeInterpolationEquallySpaced(nodes_number):
    nodes = equallySpacedNodesFromMinusPiTo2Pi(nodes_number)
    (x, y, xplt, yplt, max_err, sqr_err) = lagrangeInterpolation(nodes)

    plot.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Lagrange`a używająca ' + str(nodes_number) + ' równoodległych węzłów')
    plot.legend(('Węzły', 'wielomian Lagrange`a', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/lagrange_eqspaced_plot_%d' % nodes_number)
    plot.close()

    report.write("Lagrange %d rownoodleglych wezlow:" % nodes_number + "    max_error =%10.3E" % max_err +
                 " sqr_err =%10.3E\n" % sqr_err)

def hermiteLagrangeInterpolationChebyshev(nodes_number):
    nodes = chebyshevNodesFromMinusPiTo2Pi(nodes_number)
    (x, y, xplt, yplt) = hermiteLagrangeInterpolation(nodes)

    plot.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Hermite Lagrange używająca ' + str(nodes_number) + ' węzłów Czebyszewa')
    plot.legend(('Węzły', 'wielomian Hermite Lagrange', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/hermite_lagrange_chebyshev_plot_%d' % nodes_number)
    plot.close()

    report.write("Hermite Lagrange %d wezlow Czebyszewa:" % nodes_number)
    #+ "         max_error =%10.3E" % max_err +" sqr_err =%10.3E\n" % sqr_err)


def hermiteLagrangeInterpolationEquallySpaced(nodes_number):
    nodes = equallySpacedNodesFromMinusPiTo2Pi(nodes_number)
    (x, y, xplt, yplt) = hermiteLagrangeInterpolation(nodes)

    plot.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Hermite Lagrange używająca ' + str(nodes_number) + ' równoodległych węzłów')
    plot.legend(('Węzły', 'wielomian Hermite Lagrange', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/hermite_lagrange_eqspaced_plot_%d' % nodes_number)
    plot.close()

    report.write("Hermite Lagrange %d rownoodleglych wezlow:" % nodes_number)
    #+ "    max_error =%10.3E" % max_err + " sqr_err =%10.3E\n" % sqr_err)


def newtonCoefficients(x, y):
    n = len(y)
    arr = numpy.zeros([n, n])
    arr[:, 0] = y
    coefficients = [y[0]]

    for j in range(1, n):
        for i in range(j, n):
            arr[i][j] = (arr[i][j - 1] - arr[i - 1][j - 1]) / (x[i] - x[i - j])
            if i == j:
                coefficients.append(arr[i][j])

    return coefficients

def hermiteCoefficients(x, y):
    n = len(y)
    arr = numpy.zeros([2 * n, 2 * n])
    coefficients = [y[0]]

    for i in range(n):
        arr[2 * i, 0] = y[i]
        arr[2 * i + 1, 0] = y[i]
        arr[2 * i + 1, 1] = valueOfDerivativeInPoint(x[2 * i])
        if i > 0:
            arr[2 * i, 1] = (arr[2 * i][0] - arr[2 * i - 1][0]) / (x[2 * i] - x[2 * i - 1])

    coefficients.append(arr[1][1])

    for j in range(2, 2 * n):
        for i in range(j, 2 * n):
            arr[i][j] = (arr[i][j - 1] - arr[i - 1][j - 1]) / (x[i] - x[i - j])
            if i == j:
                coefficients.append(arr[i][j])

    return coefficients


def newtonPolynomialValue(coefficients, x, x_plt):
    n = len(x)
    y_polyn = []
    max_err = 0
    for x_polyn in x_plt:
        y_fun = valueOfFunctionInPoint(x_polyn)
        val = coefficients[0]
        for k in range(1, n):
            to_add = coefficients[k]
            for j in range(k):
                to_add *= x_polyn - x[j]
            val += to_add
        y_polyn.append(val)
        max_err = max(abs(y_fun - val), max_err)

    return y_polyn, max_err


def hermitePolynomialValue(coefficients, x, x_plt):
    n = len(x)
    y_polyn = []
    max_err = 0
    for x_polyn in x_plt:
        y_fun = valueOfFunctionInPoint(x_polyn)
        val = coefficients[0]
        for k in range(1, n):
            to_add = coefficients[k]
            for j in range(k):
                to_add *= x_polyn - x[j]
            val += to_add
        y_polyn.append(val)
        max_err = max(abs(y_fun - val), max_err)

    return y_polyn, max_err


def newtonInterpolationChebyshev(nodes_number):
    x = chebyshevNodesFromMinusPiTo2Pi(nodes_number)
    y = valuesOfFunctionInNodes(x)
    xplt = numpy.linspace(-math.pi, 2 * math.pi, num=900)
    (yplt, max_err) = newtonPolynomialValue(newtonCoefficients(x, y), x, xplt)

    y_fun = valuesOfFunctionInNodes(xplt)
    sqr_err = numpy.square(numpy.subtract(y_fun, yplt)).mean()
    sqr_err = math.sqrt(sqr_err)

    plot.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Newtona używająca ' + str(nodes_number) + ' węzłów Czebyszewa')
    plot.legend(('Węzły', 'wielomian Newtona', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/newton_chebyshev_plot_%d' % nodes_number)
    plot.close()

    report.write("Newton %d wezlow Czebyszewa:" % nodes_number + "           max_error =%10.3E" % max_err +
                 " sqr_err =%10.3E\n" % sqr_err)


def newtonInterpolationEquallySpaced(nodes_number):
    x = equallySpacedNodesFromMinusPiTo2Pi(nodes_number)
    y = valuesOfFunctionInNodes(x)
    xplt = numpy.linspace(-math.pi, 2 * math.pi, num=900)
    (yplt, max_err) = newtonPolynomialValue(newtonCoefficients(x, y), x, xplt)

    y_fun = valuesOfFunctionInNodes(xplt)
    sqr_err = numpy.square(numpy.subtract(y_fun, yplt)).mean()
    sqr_err = math.sqrt(sqr_err)

    plot.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Newtona używająca ' + str(nodes_number) + ' równoodległych węzłów')
    plot.legend(('Węzły', 'wielomian Newtona', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/newton_eqspaced_plot_%d' % nodes_number)
    plot.close()

    report.write("Newton %d rownoodleglych wezlow:" % nodes_number + "     max_error =%10.3E" % max_err +
                 " sqr_err =%10.3E\n" % sqr_err)

def hermiteInterpolationChebyshev(nodes_number):
    arr = chebyshevNodesFromMinusPiTo2Pi(nodes_number)
    y = valuesOfFunctionInNodes(arr)
    x = [0] * (nodes_number * 2)
    for i in range(nodes_number):
        x[2 * i] = arr[i]
        x[2 * i + 1] = arr[i]
    xplt = numpy.linspace(-math.pi, 2 * math.pi, num=900)
    (yplt, max_err) = hermitePolynomialValue(hermiteCoefficients(x, y), x, xplt)

    y_fun = valuesOfFunctionInNodes(xplt)
    sqr_err = numpy.square(numpy.subtract(y_fun, yplt)).mean()
    sqr_err = math.sqrt(sqr_err)

    plot.plot(arr, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Hermite używająca ' + str(nodes_number) + ' węzłów Czebyszewa')
    plot.legend(('Węzły', 'wielomian Hermite', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/hermite_chebyshev_plot_%d' % nodes_number)
    plot.close()

    report.write("Hermite %d wezlow Czebyszewa:" % nodes_number + "           max_error =%10.3E" % max_err +
                 " sqr_err =%10.3E\n" % sqr_err)


def hermiteInterpolationEquallySpaced(nodes_number):
    arr = equallySpacedNodesFromMinusPiTo2Pi(nodes_number)
    y = valuesOfFunctionInNodes(arr)
    x = [0] * (nodes_number * 2)
    for i in range(nodes_number):
        x[2 * i] = arr[i]
        x[2 * i + 1] = arr[i]
    xplt = numpy.linspace(-math.pi, 2 * math.pi, num=900)
    (yplt, max_err) = hermitePolynomialValue(hermiteCoefficients(x, y), x, xplt)

    y_fun = valuesOfFunctionInNodes(xplt)
    sqr_err = numpy.square(numpy.subtract(y_fun, yplt)).mean()
    sqr_err = math.sqrt(sqr_err)

    plot.plot(arr, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plot.xlabel('x')
    plot.ylabel('y')
    plot.title('Interpolacja Hermite używająca ' + str(nodes_number) + ' równoodległych węzłów')
    plot.legend(('Węzły', 'wielomian Hermite', 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
    plot.ylim(-4, 9)
    plot.savefig('plots/hermite_eqspaced_plot_%d' % nodes_number)
    plot.close()

    report.write("Hermite %d rownoodleglych wezlow:" % nodes_number + "     max_error =%10.3E" % max_err +
                 " sqr_err =%10.3E\n" % sqr_err)


if __name__ == '__main__':
    report = open("wyniki_interpolacji.txt", "a")

    for arg in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
        num_of_nodes = arg
        report.write("=================== LICZBA WEZLOW: %d\n" % num_of_nodes)
        #lagrangeInterpolationChebyshev(num_of_nodes)
        #lagrangeInterpolationEquallySpaced(num_of_nodes)
        #newtonInterpolationChebyshev(num_of_nodes)
        #newtonInterpolationEquallySpaced(num_of_nodes)
        #hermiteInterpolationChebyshev(num_of_nodes)
        #hermiteInterpolationEquallySpaced(num_of_nodes)
        hermiteLagrangeInterpolationChebyshev(num_of_nodes)
        hermiteLagrangeInterpolationEquallySpaced(num_of_nodes)
        report.write('\n')
