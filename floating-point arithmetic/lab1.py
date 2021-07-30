import math
import sys


def formula_1(x, y):
    return math.sqrt(x ** 2 + y ** 2)


def formula_2(v, w):
    return v * math.sqrt(1 + (w / v) ** 2)


def formula_3(v, w):
    return 2 * v * math.sqrt(1 / 4 + (w / (2 * v)) ** 2)


x = float(sys.argv[1])
y = float(sys.argv[2])
print(x, y)

v = max(abs(x), abs(y))
w = min(abs(x), abs(y))

print(formula_1(x, y))
print(formula_2(v, w))
print(formula_3(v, w))

