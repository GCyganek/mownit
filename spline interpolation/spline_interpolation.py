import sys
import math
import numpy
import matplotlib.pyplot as plt


def createPlotAndWriteToReport(xplt, yplt, x, y, plot_label, plot_file_name, interpolation_polyn_name, report_label):
    plt.plot(x, y, 'ro', xplt, yplt, 'b-', xplt, valuesOfFunctionInNodes(xplt), 'm-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(plot_label)
    plt.legend(('Węzły', interpolation_polyn_name, 'zadana funkcja f(x)'), loc='lower center', fontsize='x-small')
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


def quadraticSplineDerivativesAtNodes(z, nodes, nodes_values, from_the_left):
    if from_the_left:
        for x in range(1, len(nodes)):
            z[x] = -z[x - 1] + 2 * (nodes_values[x] - nodes_values[x - 1]) / (nodes[x] - nodes[x - 1])
    else:
        for x in range(len(nodes) - 2, -1, -1):
            z[x] = -z[x + 1] + 2 * (nodes_values[x + 1] - nodes_values[x]) / (nodes[x + 1] - nodes[x])


def quadraticSplineInterpolation(num_of_nodes, is_condition_on_left, is_linear):
    nodes = equallySpacedNodesFromMinusPiTo2Pi(num_of_nodes)
    nodes_values = valuesOfFunctionInNodes(nodes)

    z = [0] * num_of_nodes

    if is_linear:
        if is_condition_on_left:
            z[0] = (nodes_values[1] - nodes_values[0]) / (nodes[1] - nodes[0])
        else:
            z[-1] = (nodes_values[-1] - nodes_values[-2]) / \
                    (nodes[-1] - nodes[-2])
    else:
        if is_condition_on_left:
            z[0] = valueOfDerivativeInPoint(-math.pi)
        else:
            z[num_of_nodes - 1] = valueOfDerivativeInPoint(2 * math.pi)

    quadraticSplineDerivativesAtNodes(z, nodes, nodes_values, is_condition_on_left)

    xplt = equallySpacedNodesFromMinusPiTo2Pi(900)
    yplt = []
    z_index = 0

    function_values = []
    max_err = 0

    for x in xplt:
        if x > nodes[z_index + 1]:
            z_index += 1

        y = ((z[z_index + 1] - z[z_index]) / (2 * (nodes[z_index + 1] - nodes[z_index]))) * \
            math.pow(x - nodes[z_index], 2) + z[z_index] * (x - nodes[z_index]) + nodes_values[z_index]

        function_value = valueOfFunctionInPoint(x)
        max_err = max(abs(function_value - y), max_err)
        yplt.append(y)
        function_values.append(function_value)

    sqr_err = numpy.square(numpy.subtract(function_values, yplt)).mean()
    sqr_err = math.sqrt(sqr_err)

    polynom_name = "F. sklejane 2-go stopnia"
    if is_linear:
        if is_condition_on_left:
            plot_name = "Interpolacja f. sklejanymi 2-go st. - Q0''(x0) = 0 - %d węzłów" % num_of_nodes
            file_name = "plots/2stopnia/q''0/quadratic_spline_q''x0_%d" % num_of_nodes
            to_report = "Funkcje sklejane 2-go stopnia Q_0''(x0) = 0 %d węzłów:" % num_of_nodes + \
                        "         max_error =%10.3E" % max_err + " sqr_err =%10.3E\n" % sqr_err
        else:
            plot_name = "Interpolacja f. sklejanymi 2-go st. - Q(n-1)''(xn) = 0 - %d węzłów" % num_of_nodes
            file_name = "plots/2stopnia/q''n/quadratic_spline_q''x0_%d" % num_of_nodes
            to_report = "Funkcje sklejane 2-go stopnia Q_(n-1)''(xn) = 0 %d węzłów:" % num_of_nodes + \
                        "         max_error =%10.3E" % max_err + " sqr_err =%10.3E\n" % sqr_err
    else:
        if is_condition_on_left:
            plot_name = "Interpolacja f. sklejanymi 2-go st. - Q0'(x0) = f'(x0) - %d węzłów" % num_of_nodes
            file_name = "plots/2stopnia/q'0/quadratic_spline_q'x0_%d" % num_of_nodes
            to_report = "Funkcje sklejane 2-go stopnia Q_0'(x0) = f'(x0) %d węzłów:" % num_of_nodes + \
                        "         max_error =%10.3E" % max_err + " sqr_err =%10.3E\n" % sqr_err
        else:
            plot_name = "Interpolacja f. sklejanymi 2-go st. - Q(n-1)'(xn) = f'(xn) - %d węzłów" % num_of_nodes
            file_name = "plots/2stopnia/q'n/quadratic_spline_q'xn_%d" % num_of_nodes
            to_report = "Funkcje sklejane 2-go stopnia Q_(n-1)'(xn) = f'(xn) %d węzłów:" % num_of_nodes + \
                        "         max_error =%10.3E" % max_err + " sqr_err =%10.3E\n" % sqr_err

    createPlotAndWriteToReport(xplt, yplt, nodes, nodes_values, plot_name, file_name, polynom_name, to_report)


# Qk(x) = ak + bk * (x - xk) + ck * (x - xk) ^ 2 + dk * (x - xk) ^ 3 for xk <= x <= x(k + 1)
def cubicSplineInterpolation(num_of_nodes, is_condition_with_derivative):
    nodes = equallySpacedNodesFromMinusPiTo2Pi(num_of_nodes)
    a = valuesOfFunctionInNodes(nodes)

    h = nodes[1] - nodes[0]

    alpha = [0] * num_of_nodes

    l = [1] * num_of_nodes
    mi = [0] * num_of_nodes
    z = [0] * num_of_nodes

    c = [0] * num_of_nodes
    b = [0] * num_of_nodes
    d = [0] * num_of_nodes

    if is_condition_with_derivative:
        first_deriv_at_x0 = valueOfDerivativeInPoint(nodes[0])
        first_deriv_at_xn = valueOfDerivativeInPoint(nodes[-1])

        alpha[0] = 3 * (a[1] - a[0]) / h - 3 * first_deriv_at_x0
        alpha[-1] = 3 * first_deriv_at_xn - 3 * (a[-1] - a[-2]) / h

        l[0] = 2 * h
        mi[0] = 0.5
        z[0] = alpha[0] / l[0]

    for i in range(1, num_of_nodes - 1):
        alpha[i] = (3 / h) * (a[i + 1] - 2 * a[i] + a[i - 1])

    for i in range(1, num_of_nodes - 1):
        l[i] = 2 * (nodes[i + 1] - nodes[i - 1]) - h * mi[i - 1]
        mi[i] = h / l[i]
        z[i] = (alpha[i] - h * z[i - 1]) / l[i]

    if is_condition_with_derivative:
        l[-1] = h * (2 - mi[-2])
        z[-1] = (alpha[-1] - h * z[-2]) / l[-1]
        c[-1] = z[-1]

    for j in range(num_of_nodes - 2, -1, -1):
        c[j] = z[j] - mi[j] * c[j + 1]
        b[j] = (a[j + 1] - a[j]) / h - h * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h)

    xplt = equallySpacedNodesFromMinusPiTo2Pi(900)
    yplt = []
    index = 0

    function_values = []
    max_err = 0

    for x in xplt:
        if x > nodes[index + 1]:
            index += 1

        y = a[index] + b[index] * (x - nodes[index]) + c[index] * math.pow(x - nodes[index], 2) + \
            d[index] * math.pow(x - nodes[index], 3)

        function_value = valueOfFunctionInPoint(x)
        max_err = max(abs(function_value - y), max_err)
        yplt.append(y)
        function_values.append(function_value)

    sqr_err = numpy.square(numpy.subtract(function_values, yplt)).mean()
    sqr_err = math.sqrt(sqr_err)

    if is_condition_with_derivative:
        plot_name = "Interpolacja f. sklejanymi 3-go st. - Q0'(x0) = f'(x0),\nQ(n-1)'(xn) = f'(xn) - %d węzłów" % num_of_nodes
        file_name = "plots/3stopnia/q'0 q'n/cubic_spline_q'x0_q'xn_%d" % num_of_nodes
        to_report = "Funkcje sklejane 3-go stopnia Q0'(x0) = f'(x0), Q(n-1)'(xn) = f'(xn) - %d węzłów:" % num_of_nodes + \
                    "         max_error =%10.3E" % max_err + " sqr_err =%10.3E\n" % sqr_err
        polynom_name = "F. sklejane 3-go stopnia"

    else:
        plot_name = "Interpolacja f. sklejanymi 3-go st. - Q(n-1)''(xn) = Q0''(x0) = 0\n- %d wezlow" % num_of_nodes
        file_name = "plots/3stopnia/q''0 q''n/cubic_spline_q''x0_q''xn_%d" % num_of_nodes
        to_report = "Funkcje sklejane 3-go stopnia Q(n-1)''(xn) = Q0''(x0) = 0 %d wezlow:" % num_of_nodes + \
                    "         max_error =%10.3E" % max_err + " sqr_err =%10.3E\n" % sqr_err
        polynom_name = "F. sklejane 3-go stopnia"

    createPlotAndWriteToReport(xplt, yplt, nodes, a, plot_name, file_name, polynom_name, to_report)


if __name__ == '__main__':
    report = open("wyniki_interpolacji.txt", "a")

    for arg in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
        num_of_nodes = arg
        report.write("=================== LICZBA WEZLOW: %d\n" % num_of_nodes)
        quadraticSplineInterpolation(num_of_nodes, True, True)
        quadraticSplineInterpolation(num_of_nodes, False, True)
        quadraticSplineInterpolation(num_of_nodes, True, False)
        quadraticSplineInterpolation(num_of_nodes, False, False)
        cubicSplineInterpolation(num_of_nodes, False)
        cubicSplineInterpolation(num_of_nodes, True)
        report.write('\n')
