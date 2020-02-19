# Calculate cloud conduction in a hot ambient medium
# Reference: Cowie & McKee (1977)

# Everything is carried out in c.g.s unit system

from astroconst import pc, ac
from numpy import exp, log, pi, sqrt
from scipy.optimize import bisect

print "Compiling ... conduction.py ... DONE"

# PHI_S = 1.1 # Fully ionized gas, Ti = Te
# M_S = 1.5106 # Ms(1 + Ms^2/5) = 2 Phi_s; Ms ~ 2 Phi_s << 1

PHI_S = 1.0 # Fully ionized gas, Ti = Te
M_S = 1.4234 # Ms(1 + Ms^2/5) = 2 Phi_s; Ms ~ 2 Phi_s << 1
GAMMA = 5./3.
# Hint: This could be tabulated if we want different MS

class parameters():
    def __init__(self, m_s):
        self.m_s = m_s
        self.mssq = m_s ** 2
        self.phi_s = m_s * (1.0 + self.mssq / 5.0) / 2.0
        h = (11.0 + self.mssq) / (1.0 + self.mssq)
        self.ys1 = h / (h - 1)
        # H = h ** h / (h-1) ** (h-1)
        self.H = self.ys1 ** h / (self.ys1 - 1.0)
        if(self.H < 11.5): self.H = 11.5 # Approximate for supersonic case
        self.fsigma_fac = self.H ** (1.0 + self.mssq) * exp(-2.5 * self.mssq)
        self.fsigma_fac = 2.0 * self.fsigma_fac ** (1./(6. + self.mssq))
        self.sigma_onset = (self.fsigma_fac / 2.0) ** (1.2 + self.mssq / 5.0)
        self.p_fac = self.H * exp(1.8 + 0.8 * self.mssq)
        self.p_fac = self.p_fac ** (self.mssq / (6.0 + self.mssq))
        self.p_fac *= self.ys1 ** (-2.0 * self.mssq / (1.0 + self.mssq))        
        
        print "Parameter List: "
        print "--------------------------------"
        print "phi_s = ", self.phi_s
        print "ys1 = ", self.ys1
        print "h, H = ", h, self.H
        print "Fsigma_fac = ", self.fsigma_fac
        print "sigma_onset = ", self.sigma_onset
        print "p_fac = ", self.p_fac

def init_paramters(): # Hopefully I will add this later
    return 0

params = parameters(M_S)

def coulomb_log(n_e, T_e):
    if(T_e > 4.2e5):
        return 15.88 + log(T_e / sqrt(n_e))
    # 15.88 = 29.7 - log(1.e6)    
    else: # What if T_e < 4.2e5?
        return 15.88 + log(T_e / sqrt(n_e))        

# CM77EQ1; Spitzer 1962    
def kappa_cond(n_e, T_e):
    y = 1.84e-5 * T_e ** (2.5)
    y = y / coulomb_log(n_e, T_e)
    return y

# CM77EQ5
def mfp(n_e, T_e):
    y = 1.31 * n_e * pc.k * sqrt(pc.k * T_e / pc.me)
    return kappa_cond(n_e, T_e) / y

def q_sat(n_e, T_e):
    qsat = 0.34 * n_e * pc.k * T_e * sqrt(pc.k * T_e / pc.me)
    return qsat

# CM77, Between EQ31 and EQ32
# The saturated conduction parameter (sigma_0)
def sigma_cond(r_c, n_e, T_e):
    y = 1.84 / params.phi_s
    y *= mfp(n_e, T_e) / r_c
    return y

# CM77EQ54, EQ56, EQ48
# Don't use it. It works NOT well.
def sigma_onset():
    y = exp(-0.5 * params.mssq)
    y *= params.ys1 ** (2.2) / (params.ys1 - 1.0) ** (0.2)
    print "Phi_s, Ms = ", params.phi_s, params.m_s
    print "h, ys1, sigma_onset = ", h, params.ys1, y
    return y

# CM77EQ61, EQ62
def Fsigma(sigma0):
    if(sigma0 < params.sigma_onset): # Classical; approx. sigma_onset = 1.0
        return 2.0 * sigma0
    else:
        power_idx = (1.0 + params.mssq) / (6.0 + params.mssq)
        return params.fsigma_fac * sigma0 ** power_idx

# CM77EQ60
# Mass loss rate: 4*pi*Rc^2*rho_f*c_f*phis*Fsigma
def mass_loss_rate(r_c, n_h, T_h):
    n_e = n_h
    rho_h = n_h * (0.62 * pc.mh) # Completely ionized
    cs_h = sqrt(GAMMA * pc.k * T_h / (0.62 * pc.mh)) # Sound speed Right?
    y = 4.0 * pi * params.phi_s
    y *= rho_h * cs_h * r_c * r_c
    y *= Fsigma(sigma_cond(r_c, n_e, T_h))
    return y

# The evaporation timescale
def tau_evap(M_c, r_c, n_h, T_h):
    return M_c / mass_loss_rate(r_c, n_h, T_h)

def tau_evap_bs16(chi, mach, A = 0.01):
    # Be warned, this ignores the g component
    g, gm1, gp1 = 5./3., 2./3.,  8./3.
    msq = mach ** 2
    x = (gm1 * msq + 2) * (2. * g * msq - gm1) / (4. * gp1 ** 2 * msq)
    print "x = ", x
    if(x < 1): x = 1
    return 1./(A * x * sqrt(chi))

# CM77EQ22
def tau_evap_classical(r_c, n_c, n_h, T_h):
    y = 3.3e20 * coulomb_log(n_h, T_h) / 30.0
    y *= (r_c / ac.kpc * 1.e3) ** 2 * n_c * T_h ** (-2.5)
    y *= ac.yr
    return y

def tau_cloud_crushing(r_c, n_c, n_h, vrel):
    return sqrt(n_c / n_h) * r_c / vrel

def evap_pressure_fac(r_c, n_e, T_h):
    sigma0 = sigma_cond(r_c, n_e, T_h)
    # y = 3.59 * sigma0 ** 0.4
    if(sigma0 >= 1.6): y = 5.82 * (sigma0/1.6) ** 0.4
    else: y = 5.82 * (sigma0/1.6)
    if(y < 1): y = 1
    return y

def field_length(n_c, n_e, T_e, Lambda_Cool):
    y = kappa_cond(n_e, T_e)
    return sqrt(y * T_e / Lambda_Cool) / n_c

def eqn_for_q(x, *args):
    mach, phis = args[0], args[1]
    b = 1./(GAMMA * mach * mach)
    A = (8. * x - 5. * (1. + b)) ** 2
    q = 10. * phis * (1. + b - x) ** (1.5) * sqrt(x)
    y = A - (9. + 16. * q + 5. * b * (5. * b - 6.))
    return y

def solve_for_q(mach, phis=PHI_S):
    x = bisect(eqn_for_q, 0.001, 0.25, args=(mach, phis))
    b = 1./(GAMMA * mach * mach)    
    q = 10. * phis * (1. + b - x) ** (1.5) * sqrt(x)
    print x
    return q

# WARNING: This is defined relative to RH condition at M ~ infinity!
def fac_conduction_density(q, mach=0): # This is the right answer
    if(mach == 0): beta = 0.0
    else: beta = 1./(GAMMA * mach * mach)
    y = sqrt(9.+16.*q+5.*beta*(5.*beta-6.))
    return (5.*(1.+beta) + y) / (8.*(1.-q+5.*beta))

def fac_conduction_density_2(q, mach=0): # I don't remember why it's here, but it's wrong
    if(mach == 0): beta = 0.0
    else: beta = 1./(GAMMA * mach * mach)
    y = sqrt(9.+16.*q+5.*beta*(5.*beta+10.))
    return (5.*(1.+beta) + y) / (8.*(1.-q))

def fac_conduction_x(q, mach=0): 
    if(mach == 0): beta = 0.0
    else: beta = 1./(GAMMA * mach * mach)
    y = sqrt(9.+16.*q+5.*beta*(5.*beta-6.))
    return (5.*(1.+beta) - y) / 8.0

def fac_conduction_x_2(q, mach=0): # Simplified, but wrong
    if(mach == 0): beta = 0.0
    else: beta = 1./(GAMMA * mach * mach)
    y = sqrt(9.+16.*q+5.*beta*(5.*beta+10.))
    return (5.*(1.+beta) - y) / 8.0

# WARNING: This is defined relative to RH condition at M ~ infinity!
def fac_conduction_temperature(q, mach):
    if(mach == 0): beta = 0.0
    else: beta = 1./(GAMMA * mach * mach)
    y = sqrt(9.+16.*q+5.*beta*(5.*beta-6.))
    return 0.5 - 4.*q/3. + (1.+beta)*y/6. + 5.*beta*(1.+beta/6.)

def ratio_temperature(q, mach):
    x = fac_conduction_x(q, mach)
    if(mach == 0): beta = 0.0
    else: beta = 1./(GAMMA * mach * mach)
    return (1. + beta - x) / beta * x

def ratio_density(q, mach):
    x = fac_conduction_x(q, mach)
    return 1. / x

def ratio_pressure(q, mach):
    x = fac_conduction_x(q, mach)
    return ratio_temperature(q, mach) / x

def fac_conduction_pressure(q, mach):
    return fac_conduction_density(q, mach) * fac_conduction_temperature(q, mach)

def show_heat_flux(mach, M_c, v, n_a, T_a, qguess):
    import bowshock
    # conduction.show_heat_flux(3.8, 6.7e4*ac.msolar, 1.e8, 1./300., 3.e6, 0.95)
    n_e = fac_conduction_density(qguess, mach) * n_a * bowshock.rh_density_ratio(mach)
    T_e = fac_conduction_temperature(qguess, mach) * T_a * bowshock.rh_temperature_ratio(mach)
    n_c = n_e * T_e / 1.e4
    r_c = (M_c / (pi * n_c * 0.62 * pc.mh)) ** (1./3.)
    qclass = (5.6e-7 * T_e ** 2.5) * (T_e - T_a) / r_c
    qsat = 0.34 * n_e * pc.k * T_e * sqrt(pc.k * T_e / pc.me)
    qflow = 0.5 * n_a * v ** 3 * pc.mh * 0.60
    qratio = min(qsat, qclass)/qflow
    print "Density, Temperature Ratios: %5.3f, %5.3f" % (n_e/n_a, T_e/T_a)
    print "Cloud Radius: %5.3f [pc]" % (r_c * 1.e3 / ac.kpc)
    print "Classical Flux: %5.3e" % (qclass)
    print "Saturated Flux: %5.3e" % (qsat)
    print "Flow Heat Flux: %5.3e" % (qflow)
    print "Min(qcond/qflow) = %5.3f" % (qratio)
    print "Density Enhancement = %5.3f" % (fac_conduction_density(qratio, mach))
    print "Temperature Enhancement = %5.3f" % (fac_conduction_temperature(qratio, mach))    

def sanity_check():
    T_f = 3.0e6 # K
    n_f = 0.001 # cm^-3
    r_c = 0.01 * ac.kpc # SENSITIVE
    M_c = 1.e4 * ac.msolar
    rho_c = M_c / (4.18879 * r_c ** 3)
    n_c = rho_c / (pc.mh * 1.30)
    print "n_c = ", n_c
    print "t_evap = ", tau_evap(M_c, r_c, n_f, T_f) / (1.e6 * ac.yr), "Myr"
    print "t_evap (EQ64) = ", 2.8 * (n_c / n_f) * (r_c * 1.e3/ ac.kpc) / sqrt(T_f * params.phi_s * Fsigma(sigma_cond(r_c, n_f, T_f))), "Myr"
    print "t_evap_class = ", tau_evap_classical(r_c, n_c, n_f, T_f)/ac.yr/1.e6, "Myr"
    print "mfp = ", mfp(n_f, T_f) / (ac.kpc), "kpc"
    print "sigma_cond = ", sigma_cond(r_c, n_f, T_f)

def plotting():
    import matplotlib.pyplot as plt
    from numpy import linspace, logspace, log10
    # x: r_c, T_f, n_f
    # y: sigma_cond, tau_evap, tau_evap_class
    T_f = 3.0e6 # K
    n_f = 0.001 # cm^-3
    r_c = 0.01 * ac.kpc # SENSITIVE
    M_c = 1.e4 * ac.msolar
    rho_c = M_c / (4.18879 * r_c ** 3)
    n_c = rho_c / (pc.mh * 1.30)
    y1, y2, y3, y4 = [], [], [], []
    # x = logspace(6.0, 7.0, 20)
    # for x in T_f:
    #     y1.append(sigma_cond(r_c, n_f, x))
    #     y2.append(tau_evap(M_c, r_c, n_f, x) / (1.e6 * ac.yr))
    #     y3.append(tau_evap_classical(r_c, n_c, n_f, x) / (1.e6 * ac.yr))
    x = linspace(0.01*ac.kpc, 0.1*ac.kpc, 20)
    for r_c in x:
        rho_c = M_c / (4.18879 * r_c ** 3)
        n_c = rho_c / (pc.mh * 1.30)
        y1.append(sigma_cond(r_c, n_f, T_f))
        y2.append(tau_evap(M_c, r_c, n_f, T_f) / (1.e6 * ac.yr))
        y3.append(tau_evap_classical(r_c, n_c, n_f, T_f) / (1.e6 * ac.yr))
        y4.append(mfp(n_f, T_f) / (ac.kpc / 1.e3))
    plt.plot(x, y1, "b-")
    plt.plot(x, y2, "r-")
    plt.plot(x, y3, "g-")
    plt.plot(x, y4, "k-")    
    # plt.xscale("log")
    plt.yscale("log")    
    plt.show()

def plot_fac_conduction():
    import matplotlib.pyplot as plt
    q = linspace(0., 1., 100)
    y1, y2 = [], []
    for q0 in q:
        y1.append(fac_conduction_density(q0))
        y2.append(fac_conduction_temperature(q0))        
    plt.plot(q, y1, "b-")
    plt.plot(q, y2, "r-")
    plt.show()

def evaluate_saturation_point(T_star, n_ps, T_ps, L_T, T_c = 1.e4, f_R = 1.0):
    if(T_star == 0): T_star = T_ps
    a = 2.4e4 * (T_ps ** 2.5 - T_star ** 2.5) * (T_c / T_star) ** 1.53
    b = n_ps * T_ps * T_star ** (-0.5) * L_T * f_R
    return a / b - 1.0


def solve_for_saturation_point(n_ps, T_ps, L_T, T_c=1.e4, verbose=True):
    if(evaluate_saturation_point(T_c, n_ps, T_ps, L_T, T_c=T_c) <= 0):
        return -1 # No solution
    else:
        T_star = bisect(evaluate_saturation_point, T_c, T_ps, args=(n_ps, T_ps, L_T))
        if(verbose == True):
            print "classical mlr = %5.3f" % (4.46e-15 * (T_ps ** 2.5 - T_star ** 2.5))
            print "saturated mlr = %5.3f" % (1.42e-25 * n_ps * T_ps * L_T * T_star ** 1.03)
    return T_star

# Test:
def test():
    from numpy import linspace
    import matplotlib.pyplot as plt
    n_ps, L_T = 0.001, 10. * ac.pc
    T_star = 1.e4
    Tas = linspace(1.e6, 3.e7, 100)
    y = []
    for T_ps in Tas:
        # y.append(evaluate_saturation_point(T_star, n_ps, T_ps, L_T))
        y.append(solve_for_saturation_point(n_ps, T_ps, L_T))
    plt.plot(Tas, y, "b.-")
    plt.axhline(0.0)
    plt.show()
    
# Characteristic scale below which KH is suppressed
def lambda_kh(mach, chi, T_a, n_a, fcond=1.0):
    y = fcond * 23. * ac.kpc
    y *= 1.0 / (3.0 + mach * mach) / mach
    y *= sqrt(chi / 100.) * (T_a/1.e7)**2 * (0.01/n_a)
    return y

# OBSOLETE FUNCTIONS

# def evaluate_saturation_point(T_star, n_ps, T_ps, L_T):
#     a = 6.1e-7 * (T_ps ** 3.5 - T_star ** 3.5) / 3.5 / L_T
#     if(T_star == 0): T_star = T_ps
#     b = 1.7e-11 * n_ps * T_ps * T_star ** 0.5
#     # print a, b
#     return a / b - 1.0

# def solve_for_saturation_point(n_ps, T_ps, L_T, verbose=True):
#     if(evaluate_saturation_point(1.e4, n_ps, T_ps, L_T) <= 0):
#         return -1 # No solution
#     else:
#         T_star = bisect(evaluate_saturation_point, 1.e4, T_ps, args=(n_ps, T_ps, L_T))
#         if(verbose == True):
#             print "classical mlr = %5.3f" % (1.56e-15 * T_ps ** 2.5)
#             print "saturated mlr = %5.3f" % (5.39e-19 * n_ps * sqrt(T_star) * L_T / 3.5)
#     return T_star
