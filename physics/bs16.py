
# Some formula from BS16
# Reference: Bruggen and Scannapieco 2016

# Everything is carried out in c.g.s unit system

from astroconst import pc, ac
from numpy import exp, log, pi, sqrt, linspace
from scipy.optimize import bisect

GAMMA = 5./3.
GAMMAP1 = GAMMA + 1.0
GAMMAM1 = GAMMA - 1.0

def rh_temperature_ratio(mach):
    mssq = mach * mach
    y = GAMMAP1 + 2.0 * GAMMA * (mssq - 1.0)
    y *= (GAMMAP1 + GAMMAM1 * (mssq - 1.0))
    y /= GAMMAP1 ** 2 * mssq
    return y

def cs_iso(T):
    if(T > 2.e4): return sqrt(GAMMA * pc.k * T / 0.62 / pc.mh)
    else: return sqrt(GAMMA * pc.k * T / 1.30 / pc.mh)    

# Cloud crushing time
def tau_cloud_crushing(r_c, n_c, n_h, vrel):
    ksi = sqrt(n_c / n_h)
    return ksi * r_c / vrel

def gfunc(mach, Tevap=3.e6, A=0.01, Nc=3.e20, ksi=300.0, LambdaCool=1.e-22):
    g = 3.5 * (A/0.01) * (0.5/0.5) * (Nc/3.e20) * (3.e6 / Tevap) * mach * sqrt(1000./ksi)
    return (sqrt(1.0 + 4.0 * g) - 1.0) / (2.0 * g)

def fM(mach):
    y = rh_temperature_ratio(mach) / 4.0
    if(y < 1): y = 1
    return y

def tau_ratio(mach, A=0.01, ksi=300.0, Nc=3.e20):
    y = 1./gfunc(mach, A=A, ksi=ksi, Nc=Nc)
    return y / (A * fM(mach) * sqrt(ksi))
    
def velocity(mach, v_h, Tevap=3.e6, A=0.01, ksi=300.0, Nc=3.e20):
    tratio = tau_ratio(mach, A=A, ksi=ksi, Nc=Nc)
    c_evap = cs_iso(Tevap)
    A = 0.4 * v_h / sqrt(ksi) * sqrt(mach / 30. * tratio)
    B = 0.6 * c_evap * sqrt(mach / 30. / tratio)
    return A + B

mach = 11.4
v_h = 3000.e5
ksi = 300.0
Nc = 3.e21
print tau_ratio(mach, ksi=ksi, Nc=Nc), velocity(mach, v_h, ksi=ksi, Nc=Nc)/1.e5

def EQ33(x, mach, phis=1.1):
    A = 1.25 + 1.25 / mach**2
    B = 2.5 * phis * (1.0 - x + 1./mach**2) ** (1.5)
    C = 0.25# + 1.25 / mach**2
    return x * x - A * x - B * sqrt(x) + C

def EQMy(x, mach=1.0, phis=1.1):
#def EQMy(x, *args):
    # mach = args[0]
    # phis = args[1]
    b = 1./(GAMMA * mach * mach)
    A = (8. * x - 5. * (1. + b)) ** 2
    q = 10. * phis * (1. + b - x) ** (1.5) * sqrt(x)
    y = A - (9. + 16. * q + 5. * b * (5. * b - 6.))
    return y

def EQMy_Classical(x, mach=1.0, phis=1.1):
    b = 1./(GAMMA * mach * mach)
    tratio = 1.0
    return y

# def EQMy(x, mach, phis=1.1):
#     b = 1./(mach * mach)
#     A = (8. * x - 5. * (1. + b)) ** 2
#     q = 10. * phis * (1. + b - x) ** (1.5) * sqrt(x)
#     y = A - (9. + 16. * q + 5. * b * (5. * b + 10.))
#     return y

def solve_for_q(mach, phis=1.1):
    x = bisect(EQMy, 0.001, 0.25, args=(mach, phis))
    b = 1./(GAMMA * mach * mach)    
    q = 10. * phis * (1. + b - x) ** (1.5) * sqrt(x)
    print "x = ", x
    return q

def plotEQMy(mach):
    import matplotlib.pyplot as plt
    xr = linspace(0.001, 10.0, 1000)
    yr = EQMy(xr, mach, 1.1)
    plt.plot(xr, yr, "b-")
    plt.plot([0.001, 10.0], [0.0, 0.0], "k--")
    plt.show()

def plotEQ33(mach):
    import matplotlib.pyplot as plt
    xr = linspace(0.001, 0.1, 1000)
    yr = EQ33(xr, mach)
    plt.plot(xr, yr, "b-")
    plt.plot([0.001, 0.1], [0.0, 0.0], "k--")
    plt.show()

def q_eff(mach, x, phis=1.1):
    A = 10. * phis
    y = 1. + 1./(GAMMA*mach**2) - x
    # y = 1. + 1./(mach**2) - x    
    return A * y ** (1.5) * x ** (0.5)
    
