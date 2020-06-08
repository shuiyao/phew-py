from numpy import genfromtxt, array, exp, mean, median
import ioformat
import matplotlib.pyplot as plt
from math import factorial

f = "example.wind"
h, rc = ioformat.rcol(f, [5,6])
ratio = array(rc) / array(h)
y = (ratio ** 2) * 1000.0

def p(n, y):
    return y ** n * exp(-y) / (float)(factorial(n))

P0 = p(0, y)
P1 = p(1, y)
P2 = p(2, y)

# P(x) = y^xe^-y/x!, y = mean(x)

