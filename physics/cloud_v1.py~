# Version 1.0:
# The first version of cynlindrical cloud.
import conduction
import bowshock
from astroconst import pc, ac
from numpy import exp, log, pi, sqrt, logspace, log10
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['mathtext.default'] = "regular"
mpl.rcParams['axes.labelsize'] = "large"

RES_T = 0.01 # time resolution, dt = tau_cc * RES_T
END_T = 200. # End time, tend = tau_cc * END_T
M_S = 1.5106        
FACDENS, FACTEMP = 1.0, 1.0

class parameters():
    def __init__(self, q=0.0, mach=0.0):
        self.eta_s = conduction.fac_conduction_density(q, mach)
        self.tau_s = conduction.fac_conduction_temperature(q, mach)
        self.M_c = 6.7e4 * ac.msolar
        self.v_h = 1000. * 1.e5
        self.T_eq = 1.e4
        self.T_a = 3.e6
        self.n_c = 1.0
        self.n_a = 1.0 / 300.
        self.r_0 = 3.0 * ac.kpc
        self.f_cond = 1.0
        self.T_vir = self.T_a
        self.q_s = q

class study_case():
    def __init__(self, params):
        self.M_c = params.M_c
        self.vrel = params.v_h
        self.v_h = params.v_h
        self.T_c = params.T_eq
        self.T_a = params.T_a
        self.n_c = params.n_c
        self.n_a = params.n_a
        self.eta_s = params.eta_s
        self.tau_s = params.tau_s
        self.f_cond = params.f_cond
        self.T_vir = params.T_vir
        self.r = params.r_0
        self.q_s = params.q_s
        if(params.T_vir == 0):
            self.gravity = False
        else:
            self.gravity = True
        self.v_c = 0.0
        self.a_ram = 0.0
        self.time = 0.0
        self.t75 = []
        self.t50 = []
        self.t10 = []        
    def derived_quantities(self, verbose=False):
        if(self.T_c > 8.e3): self.mu = 0.62
        else: self.mu = 1.30
        self.rho_c = self.n_c * self.mu * pc.mh
        self.rho_a = self.n_a * self.mu * pc.mh
        self.r_c = (self.M_c / (4.18879 * self.rho_c)) ** (1./3.)
        self.cs_a = bowshock.cs_iso(self.T_a)
        self.cs_c = bowshock.cs_iso(self.T_c)
        self.p_c = self.n_c * self.T_c * pc.k
        self.mach = self.vrel / self.cs_a
        self.tau_cc = bowshock.tau_cloud_crushing(self.r_c, self.n_c, self.n_a, self.vrel)
        self.N_c = self.r_c * self.n_c
        # isothermal profile
        self.V_c = sqrt(2. * pc.k * self.T_vir / (0.6 * pc.mh))
    def cloud_shock(self):
        v_shock_approx = self.r_c / self.tau_cc
        T_ps = bowshock.rh_temperature_ratio(self.mach) * self.tau_s
        if(T_ps < 1): T_ps = 1
        T_ps = T_ps * self.T_a
        n_ps = self.n_a * bowshock.rh_density_ratio(self.mach) * self.eta_s
        p_ps = n_ps * T_ps * pc.k
        p_therm_c = self.n_c * self.T_c * pc.k
        p_therm_a = self.n_a * self.T_a * pc.k
        self.v_shock = sqrt(p_ps / self.p_c) * self.cs_c
        # WARNING: Ignore the pressure difference between stagnation point and the shock front
        # if(self.mach < 1.2):
        #     p_ps *= bowshock.postshock_pressure_fac(self.mach)
        self.rho_c = self.rho_c * p_ps / p_therm_a # Isothermal shock condition
        self.r_c = (self.M_c / self.rho_c / pi) ** (1./3.) # cylindrical geometry
        self.n_c = self.rho_c / (pc.mh * self.mu) # Assuming neutral
        self.r_c = (self.M_c / self.rho_c / pi) ** (1./3.) # cylindrical geometry
        p_cond = conduction.evap_pressure_fac(self.r_c, n_ps, T_ps) * p_therm_a
        self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps) * self.f_cond
        self.tau_khi = bowshock.fac_khi(self.mach) * self.tau_cc
        self.N_c = self.M_c / (self.mu * pc.mh * pi * self.r_c * self.r_c)
        self.L = self.r_c
        # self.vrel = self.vrel - v_shock
        # fac_ram = bowshock.ram_pressure_fac(self.mach) * params.ram
        # p_ps = fac_ram * self.rho_a * self.vrel ** 2
        self.a_ram = (p_ps - p_therm_a) * pi * (self.r_c) ** 2 / self.M_c        
    def dynamical_update(self, dt):
        self.eta_s = conduction.fac_conduction_density(self.q_s, self.mach)
        self.tau_s = conduction.fac_conduction_temperature(self.q_s, self.mach)
        T_ps = bowshock.rh_temperature_ratio(self.mach) * self.tau_s
        if(T_ps < 1): T_ps = 1
        T_ps = T_ps * self.T_a
        n_ps = self.n_a * bowshock.rh_density_ratio(self.mach) * self.eta_s
        self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps) * self.f_cond
        self.tau_evap = self.tau_evap / (self.L / (2.0 * self.r_c))
        mlr = conduction.mass_loss_rate(self.r_c, n_ps, T_ps) * 0.5 # Head
        mlr += conduction.mass_loss_rate(self.r_c, FACDENS*self.n_a, FACTEMP*self.T_a) * self.L / (2.0 * self.r_c) # Tail
        self.tau_evap = self.M_c / mlr * self.f_cond
        # Calculate Acceleration
        p_therm_c = self.n_c * self.T_c * pc.k
        p_therm_a = self.n_a * self.T_a * pc.k
        p_ps = n_ps * T_ps * pc.k
        self.a_ram = (p_ps - p_therm_a) * pi * self.r_c ** 2 / self.M_c
        self.tau_khi = bowshock.fac_khi(self.mach) * self.tau_cc
        # Update Cloud Properties
        self.time = self.time + dt
        self.r = self.r + self.vrel * dt
        self.vrel = self.vrel - self.a_ram * dt
        if(self.gravity == True):
            a_grav = self.V_c ** 2 / self.r
            self.vrel = self.vrel - a_grav * dt
        self.mach = self.vrel / self.cs_a
        P_evap = 3.59 * conduction.sigma_cond(self.r_c, n_ps, T_ps)**0.28
        if(self.tau_evap > 0):
            if(P_evap > 1.0): tau_disrupt = self.tau_evap
            else: tau_disrupt = 1./(1./self.tau_evap + 1./self.tau_khi)
        else: tau_disrupt = self.tau_khi
        self.M_c = self.M_c * exp(-dt / tau_disrupt)
        print "Mc=%4.2f t(KHI,evap,cc)=%5.2f, %5.2f, %5.2f mach=%4.2f Pevap/Pth=%4.2f" % (self.M_c / params.M_c, self.tau_khi / ac.myr, self.tau_evap / ac.myr, self.tau_cc / ac.myr, self.mach, P_evap)
        # ---------------- The Expansion ----------------
        rho_eq = self.rho_a * self.T_a / self.T_c
        L_eq = self.M_c / (rho_eq * self.r_c * self.r_c * pi)
        v_exp = self.v_shock - (self.v_h - self.vrel)
        if(v_exp < 0): v_exp = 0
        self.L = self.L + v_exp * dt
        # Here using free expansion approximation: v = 2.0 * c / (gamma - 1)
        self.r_c = sqrt(self.M_c / (self.N_c * pc.mh * self.mu) / pi)
        if(self.L > L_eq): self.L = L_eq
        self.rho_c = self.M_c / (self.L * self.r_c * self.r_c * pi)
        self.n_c = self.rho_c / (pc.mh * self.mu)

params = parameters(q=0.75, mach=3.8)
def run(params):
    cloud = study_case(params)
    M_0 = params.M_c
    t75, t50, t10 = -1., -1., -1.
    cloud.derived_quantities()
    cloud.cloud_shock()
    tevap = cloud.tau_evap / cloud.tau_cc
    dt = cloud.tau_cc * RES_T

    ttot = 0.0
    # t, dv, r, R_c, L_c, n_c
    cloud.t75 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    cloud.t50 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    cloud.t10 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]    
    while ttot < END_T:
        cloud.dynamical_update(dt)
        if(cloud.M_c / M_0 < 0.75 and cloud.t75[0] == -1):
            cloud.t75 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
        if(cloud.M_c / M_0 < 0.5 and cloud.t50[0] == -1):
            cloud.t50 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
        if(cloud.M_c / M_0 < 0.10 and cloud.t10[0] == -1):
            cloud.t10 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
        if(cloud.M_c / M_0 < 0.1):
            print "STOP. Cloud is nearly disrupted. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
                (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
            break
        if(cloud.mach < 1):
            print "STOP. Flow becomes subsonic. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
                (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
            break
        ttot += RES_T
    print "End. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
        (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
    return cloud

print "Compiled."
