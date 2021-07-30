import sys
import math
import numpy
import matplotlib.pyplot as plt


def createPlotAndWriteToReport(xplt, yplt, x, y, plot_label, plot_file_name, approx_polyn_name, report_label):
    plt.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(plot_label)
    plt.legend(('Punkty dykretyzacji', approx_polyn_name, 'zadana funkcja f(x)'), loc='lower center',
               fontsize='x-small')
    plt.ylim(-4, 9)
    plt.savefig(plot_file_name)
    plt.close()

    report.write(report_label)


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


def equallySpacedNodesFromMinusPiTo2Pi(nodes_number):
    result = []
    space_between = (3 * math.pi) / (nodes_number - 1)
    for k in range(0, nodes_number):
        xk = -math.pi + space_between * k
        result.append(xk)
    return result


def sqrError(values, function_values):
    sqr_err = 0
    for index in range(len(values)):
        sqr_err += math.pow(values[index] - function_values[index], 2)
    sqr_err = sqr_err / len(values)
    sqr_err = math.sqrt(sqr_err)
    return sqr_err


def approxPolynomValuesUsingMonomials(coefficients, x):
    m = len(coefficients)
    values = []
    function_values = []
    max_err = 0
    for point in x:
        polynom_value = 0
        for index in range(m):
            polynom_value += coefficients[index] * math.pow(point, index)
        values.append(polynom_value)
        function_value = valueOfFunctionInPoint(point)
        max_err = max(abs(polynom_value - function_value), max_err)
        function_values.append(function_value)

    sqr_err = sqrError(values, function_values)

    return values, max_err, sqr_err


def approximationUsingMonomials(points_number, functions_number):
    points = equallySpacedNodesFromMinusPiTo2Pi(points_number)

    G = numpy.zeros((functions_number, functions_number))

    for x_index in range(functions_number):
        for y_index in range(functions_number):
            if x_index == 0 and y_index == 0:
                value = points_number
            else:
                value = 0
                for point in points:
                    value += math.pow(point, x_index + y_index)
            G[x_index][y_index] = value

    B = numpy.zeros((functions_number, 1))

    for index in range(functions_number):
        value = 0
        for point in points:
            value += valueOfFunctionInPoint(point) * math.pow(point, index)
        B[index][0] = value

    coefficients = numpy.linalg.solve(G, B)

    x = numpy.linspace(-math.pi, 2 * math.pi, 900)
    approx_polinom_values, max_err, sqr_err = approxPolynomValuesUsingMonomials(coefficients, x)

    createPlotAndWriteToReport(x, approx_polinom_values,
                               points, valuesOfFunctionInNodes(points),
                               "Aproksymacja sredniokwadratowa z uzyciem jednomianow i\n"
                               " %d punktow dyskretyzacji oraz %d funkcji bazowych" % (points_number, functions_number),
                               "monomials_plots/%d_fun_%d_pts_approx_monomials" % (functions_number, points_number),
                               "Funkcja aproksymujaca",
                               "\tmax_error=%10.3E\tsqr_err=%10.3E\n"
                               % (max_err, sqr_err))


def transform(x):
    return (2 / 3) * x - (math.pi / 3)


def akCoefficient(k, points):
    ak = (1 / (len(points) / 2))
    _sum = 0
    for point in points:
        transformed_point = transform(point)
        _sum += valueOfFunctionInPoint(point) * math.cos(k * transformed_point)
    ak *= _sum
    return ak


def bkCoefficient(k, points):
    bk = (1 / (len(points) / 2))
    _sum = 0
    for point in points:
        transformed_point = transform(point)
        _sum += valueOfFunctionInPoint(point) * math.sin(k * transformed_point)
    bk *= _sum
    return bk


def approximationUsingTrigonometricPolynomials(points_number, m):
    points = equallySpacedNodesFromMinusPiTo2Pi(points_number)

    ak_coefficients = [0] * m
    bk_coefficients = [0] * m

    for index in range(m):
        ak_coefficients[index] = akCoefficient(index, points)
        bk_coefficients[index] = bkCoefficient(index, points)

    x = numpy.linspace(-math.pi, 2 * math.pi, 900)

    values = []
    max_err = 0

    function_values = []

    for point in x:
        transformed_point = transform(point)
        value = ak_coefficients[0] / 2 + ak_coefficients[-1] * math.cos(m * transformed_point)
        for index in range(1, m - 1):
            value += ak_coefficients[index] * math.cos(index * transformed_point) + bk_coefficients[index] \
                     * math.sin(index * transformed_point)
        values.append(value)
        function_value = valueOfFunctionInPoint(point)
        function_values.append(function_value)
        max_err = max(abs(value - function_value), max_err)

    sqr_err = sqrError(values, function_values)

    createPlotAndWriteToReport(x, values,
                               points, valuesOfFunctionInNodes(points),
                               "Aproksymacja sredniokwadratowa z f trygonometrycznych i\n"
                               " %d punktow dyskretyzacji oraz %d funkcji bazowych" % (points_number, m * 2),
                               "trygonometric_plots/%d_fun_%d_pts_approx_trygonom" % (m * 2, points_number),
                               "Funkcja aproksymujaca",
                               "\tmax_error=%10.3E\tsqr_err=%10.3E\n"
                               % (max_err, sqr_err))


if __name__ == '__main__':
    report = open("wyniki_aproksymacji.txt", "a")
    if sys.argv[5] == "all":
        print_all = True
    elif sys.argv[5] == "%3":
        print_all = False
    else:
        sys.exit()

    if sys.argv[6] == "alg":
        report.write("Aproksymacja algebraiczna:\n")
        for arg in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
            for arg2 in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
                if (arg % 3 == 0 and arg2 % 3 == 0) or print_all is True:
                    num_of_points = arg
                    num_of_functions = arg2
                    if num_of_functions <= num_of_points:
                        report.write("\t liczba wezlow: %d, liczba f. bazowych: %d =>"
                                     % (num_of_points, num_of_functions))
                        approximationUsingMonomials(num_of_points, num_of_functions)

    elif sys.argv[6] == "tryg":
        report.write("Aproksymacja trygonometryczna:\n")
        for arg in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
            for arg2 in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
                if (arg % 3 == 0) or print_all is True:
                    num_of_points = arg
                    m = arg2
                    if 2 * m < num_of_points:
                        report.write("\t liczba wezlow: %d, liczba f. bazowych: %d =>"
                                     % (num_of_points, 2 * m))
                        approximationUsingTrigonometricPolynomials(num_of_points, m)
