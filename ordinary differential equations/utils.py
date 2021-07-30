import numpy as np
import math

x0 = np.pi / 6
xk = 3 * np.pi / 2
m = 3
k = 2

def function(x, y):
    return k * m * y * np.sin(m * x) + math.pow(k, 2) * m * np.sin(m * x) * np.cos(m * x)

def realFunction(x):
    return np.exp(-k * np.cos(m * x)) - k * np.cos(m * x) + 1

def realFunctionValues(x, k):
    y = np.zeros(k)
    for i in range(k):
        y[i] = realFunction(x[i])
    return y

def sqrError(values, real_values):
    sqr_err = 0
    for index in range(len(values)):
        sqr_err += math.pow(values[index] - real_values[index], 2)
    sqr_err = sqr_err / len(values)
    sqr_err = np.sqrt(sqr_err)
    return sqr_err
