import math
import sys


def functionValue(x):
    return (x - 1) * math.exp(-10 * x) + math.pow(x, 15)

def functionDervativeValue(x):
    return -10 * math.exp(-10 * x) * x + 11 * math.exp(-10 * x) + 15 * math.pow(x, 14)


def newtonRaphsonMethod1(x0, ro):
    step = 0
    report.write('Newton-Raphson method for x0 = %.1f, ro = %g, stop = step:\n' % (x0, ro))
    print('Newton-Raphson method for x0 = %.1f, ro = %g:' % (x0, ro))
    while True:
        x1 = x0 - functionValue(x0) / functionDervativeValue(x0)
        step += 1
        if abs(x1 - x0) < ro:
            break
        x0 = x1
    result = '\tResult: Iteration = ' + str(step) + ' | result = ' + str(x1) + '\n'
    report.write(result + '\n')
    print(result)

def newtonRaphsonMethod2(x0, ro):
    step = 0
    report.write('Newton-Raphson method for x0 = %.1f, ro = %g, stop = func:\n' % (x0, ro))
    print('Newton-Raphson method for x0 = %.1f, ro = %g:' % (x0, ro))
    while True:
        x1 = x0 - functionValue(x0) / functionDervativeValue(x0)
        step += 1
        if abs(functionValue(x1)) < ro:
            break
        x0 = x1
    result = '\tResult: Iteration = ' + str(step) + ' | result = ' + str(x1) + '\n'
    report.write(result + '\n')
    print(result)


if __name__ == '__main__':
    report = open("metoda_newtona.txt", "a")
    ro = float(sys.argv[1])
    condition = sys.argv[2]

    x0 = -1
    while x0 <= 0.7:
        if condition == 'step':
            newtonRaphsonMethod1(x0, ro)
        elif condition == 'func':
            newtonRaphsonMethod2(x0, ro)
        x0 += 0.1
