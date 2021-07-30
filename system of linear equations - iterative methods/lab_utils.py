import random
import numpy as np

def getRandom():  # from -100 to 100
    return (random.random() - 0.5) * 200

def buildMatrices(size):
    a = np.zeros((size, size))
    x = np.zeros(size)

    for i in range(size):
        x[i] = 1
        for j in range(size):
            if i == j:
                a[i][j] = 10
            else:
                a[i][j] = 1 / (abs(i - j) + 5)
    b = np.matmul(a, x)
    return a, x, b


def condition1(x0, x1, ro):
    if np.linalg.norm(x1 - x0) < ro:
        return True
    return False


def condition2(a, x1, b, ro):
    if np.linalg.norm(np.matmul(a, x1) - b) < ro:
        return True
    return False
