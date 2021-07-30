import math
import sys
import random

import numpy as np

def getRandom():  # from -50 to 50
    return (random.random() - 0.5) * 100

def equations_system(x1x2x3):
    x1, x2, x3 = x1x2x3
    return [x1 ** 2 + x2 ** 2 - x3 ** 2 - 1,
            x1 - 2 * x2 ** 3 + 2 * x3 ** 2 + 1,
            2 * x1 ** 2 + x2 - 2 * x3 ** 2 - 1]


def jacobian(x1x2x3):
    x1, x2, x3 = x1x2x3
    return [[2 * x1, 2 * x2, -2 * x3],
            [1, -6 * x2 ** 2, 4 * x3],
            [4 * x1, 1, -4 * x3]]


def newtonRaphsonMethod1(vec, epsilon):
    step = 0
    while True:
        J = np.array(jacobian(vec))
        F = np.array(equations_system(vec))

        diff = np.linalg.solve(J, -F)
        vec = vec + diff

        step += 1
        if np.linalg.norm(diff) < epsilon:
            break
        if step > 1000:
            print("Cant find result")
            report.write("\tCant find result\n")
            return

    result = "x1 = " + str(vec[0]) + ", x2 = " + str(vec[1]) + ", x3 = " + str(vec[2]) + " | iterations = " + str(step)
    report.write("\t" + result + "\n")
    print(result)


def newtonRaphsonMethod2(vec, epsilon):
    step = 0
    while True:
        J = np.array(jacobian(vec))
        F = np.array(equations_system(vec))

        diff = np.linalg.solve(J, -F)
        vec = vec + diff

        step += 1
        if np.linalg.norm(equations_system(vec)) < epsilon:
            break
        if step > 1000:
            print("Cant find result")
            report.write("\tCant find result\n")
            return

    result = "x1 = " + str(vec[0]) + ", x2 = " + str(vec[1]) + ", x3 = " + str(vec[2]) + " | iterations = " + str(step)
    report.write("\t" + result + "\n")
    print(result)


def getSign(x):
    if x > 0:
        return '+'
    elif x < 0:
        return '-'
    else:
        return '0'


if __name__ == '__main__':
    report = open("uklad_rownan.txt", "a")
    ro = float(sys.argv[1])
    condition = sys.argv[2]

    start = [getRandom(), getRandom(), getRandom()]
    initial_vector = "Initial vector: x = %10.3E, y = %10.3E, z = %10.3E" % (start[0], start[1], start[2])
    initial_vector_signs = "(" + getSign(start[0]) + ", " + getSign(start[1]) + ", " + getSign(start[2]) + ")"
    print(initial_vector + "\t" + initial_vector_signs)

    if condition == 'step':
        report.write("\nNewton method for solving nonlinear equations system, ro = %g, stop: step\n" % ro)
        report.write(initial_vector + "\t" + initial_vector_signs + "\n")
        newtonRaphsonMethod1(start, ro)
    elif condition == 'func':
        report.write("\nNewton method for solving nonlinear equations system, ro = %g, stop: func\n" % ro)
        report.write(initial_vector + "\t" + initial_vector_signs + "\n")
        newtonRaphsonMethod2(start, ro)
