from astroconst import pc, ac
from numpy import exp, log, pi, linspace
from scipy.optimize import bisect
import matplotlib.pyplot as plt

def func(y, x0=1.0, t0=0.0):
    return exp(x0 * (y - 1.0)) * (1.0 + y * t) - 1.0

# x0 = 2.0
# ts = linspace(0., 10, 100)
# ys = []
# for t in ts:
#     ys.append(bisect(func, 0.01, 1.0, args=(x0, t)))
# plt.plot(ts, ys, "b-")
# plt.show()

rho0 = 10. * pc.mh * 0.6
r_c = 20.*ac.pc
mlra = 4.46e-15 * 1.e7** 2 / (rho0 * r_c ** 2)
print log(mlra * 50. * ac.myr)

