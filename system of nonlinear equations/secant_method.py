import math
import sys


def functionValue(x):
    return (x - 1) * math.exp(-10 * x) + math.pow(x, 15)


def secantMethod1(x0, x1, ro):
    step = 0
    report.write('Secant method for x0 = %.1f, x1 = %.1f, ro = %g, stop = step:\n' % (x0, x1, ro))
    print('Secant method for x0 = %.1f, x1 = %.1f, ro = %g, stop = step:' % (x0, x1, ro))
    while abs(x1 - x0) >= ro:
        x2 = x1 - ((x1 - x0) / (functionValue(x1) - functionValue(x0))) * functionValue(x1)
        x0 = x1
        x1 = x2
        step += 1
    result = '\tResult: Iteration = ' + str(step) + ' | result = ' + str(x1) + '\n'
    report.write(result + '\n')
    print(result)


def secantMethod2(x0, x1, ro):
    step = 0
    report.write('Secant method for x0 = %.1f, x1 = %.1f, ro = %g, stop: func:\n' % (x0, x1, ro))
    print('Secant method for x0 = %.1f, x1 = %.1f, ro = %g, stop: func:' % (x0, x1, ro))
    while abs(functionValue(x1)) >= ro:
        x2 = x1 - ((x1 - x0) / (functionValue(x1) - functionValue(x0))) * functionValue(x1)
        x0 = x1
        x1 = x2
        step += 1
    result = '\tResult: Iteration = ' + str(step) + ' | result = ' + str(x1) + '\n'
    report.write(result + '\n')
    print(result)


if __name__ == '__main__':
    report = open("metoda_siecznych.txt", "a")
    ro = float(sys.argv[1])
    condition = sys.argv[2]

    x0 = 0.7
    x1 = 0.6
    while x1 >= -1:
        if condition == 'step':
            secantMethod1(x0, x1, ro)
        elif condition == 'func':
            secantMethod2(x0, x1, ro)
        x1 -= 0.1

    x0 = -1
    x1 = -0.9
    while x1 <= 0.7:
        if condition == 'step':
            secantMethod1(x0, x1, ro)
        elif condition == 'func':
            secantMethod2(x0, x1, ro)
        x1 += 0.1
