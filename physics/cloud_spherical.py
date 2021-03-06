# Version 1.0
# Spherical cloud. 

import conduction
import bowshock
import simplewave
from astroconst import pc, ac
from numpy import exp, log, pi, sqrt, logspace, log10
import matplotlib.pyplot as plt
import os.path
import sys

import matplotlib as mpl
mpl.rcParams['mathtext.default'] = "regular"
mpl.rcParams['axes.labelsize'] = "large"

RES_T = 0.01 # time resolution, dt = tau_cc * RES_T
END_T = 200. # End time, tend = tau_cc * END_T
M_S = 1.5106        
FACDENS, FACTEMP = 1.0, 1.0

ADJUST_RCLOUD_FOR_VAPORPRESSURE = True

pars_list = {"M_c":1, "v_h":1, "T_a":1, "T_eq":1, "n_c":1, "chi":1, "q_s":1, "mach":1, "f_cond":1}
class parameters():
    def __init__(self, modelname):
        # if(os.path.isfile(modelname+".param")): fpar = open(modelname+".param", "r")
        # else: sys.exit(modelname+".param not found!")
        try: fpar = open(modelname+".param", "r")
        except: print modelname+".param not found!"
        for line in fpar:
            spt = line.split()
            if(spt[1] in pars_list.keys()): pars_list[spt[1]] = float(spt[0])
            else:
                print "Error: parameter "+spt[1]+" not defined!"
        fpar.close()
        self.model = modelname
        self.q_s = pars_list['q_s']
        self.mach = pars_list['mach']
        self.eta_s = 1./conduction.fac_conduction_x(self.q_s, self.mach) # density ratio
        self.tau_s = conduction.ratio_temperature(self.q_s, self.mach) # temperature ratio
        self.M_c = pars_list['M_c'] * ac.msolar
        self.v_h = pars_list['v_h'] * 1.e5
        self.T_eq = pars_list['T_eq']
        self.T_a = pars_list['T_a']
        self.n_c = pars_list['n_c']
        self.n_a = self.n_c / pars_list['chi']
        self.r_0 = 3.0 * ac.kpc
        self.f_cond = pars_list['f_cond']
        self.T_vir = 0.0 # If non-zero, turn on gravity

class study_case():
    def __init__(self, params):
        self.model = params.model
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
        self.n_poly = 1.0 # By default, isothermal
        self.v_max_iso = 3.0 # The velocity cut to the isothermal profile
        self.pressure_factor = 1.0
        if(params.T_vir == 0):
            self.gravity = False
        else:
            self.gravity = True
            print "Gravity is ON."
        self.v_c = 0.0
        self.a_ram = 0.0
        self.time = 0.0
        self.tclock = 0.0
        self.L_prev = 0.0
        self.L = 0.0
        self.t75 = []
        self.t50 = []
        self.t10 = []        
    def derived_quantities(self, verbose=False):
        if(self.T_c > 8.e3): self.mu = 0.62
        else: self.mu = 1.30
        self.rho_c = self.n_c * self.mu * pc.mh
        self.rho_c0 = self.rho_c
        self.rho_a = self.n_a * self.mu * pc.mh
        self.r_c = (self.M_c / (4.18879 * self.rho_c)) ** (1./3.)
        self.r_c0 = self.r_c
        self.cs_a = bowshock.cs_iso(self.T_a)
        self.cs_c = bowshock.cs_iso(self.T_c)
        self.p_c = self.n_c * self.T_c * pc.k
        self.mach = self.vrel / self.cs_a
        self.tau_cc = bowshock.tau_cloud_crushing(self.r_c, self.n_c, self.n_a, self.vrel)
        # isothermal profile
        self.V_c = sqrt(2. * pc.k * self.T_vir / (0.6 * pc.mh))
    def cloud_shock(self):
        self.v_shock_approx = self.r_c / self.tau_cc # Good when conduction is not dynamically important
        T_ps = max(self.tau_s, 1.0)
        T_ps = T_ps * self.T_a
        n_ps = self.n_a * self.eta_s
        p_ps = n_ps * T_ps * pc.k
        p_therm_c = self.n_c * self.T_c * pc.k
        p_therm_a = self.n_a * self.T_a * pc.k
        self.v_shock = sqrt(p_ps / self.p_c) * self.cs_c
        self.rho_c = self.rho_c * p_ps / p_therm_a
        # ----------------------------------------------------------------
        # note 1: How do we define R_c after cloud shock?
        self.r_c = (0.75 * self.M_c / self.rho_c / pi) ** (1./3.) # # SPHERICAL cloud
        if(ADJUST_RCLOUD_FOR_VAPORPRESSURE):
            p_evap = 3.59 * conduction.sigma_cond(self.r_c, n_ps, T_ps)**0.28
            self.rho_c *= p_evap
            rc_prev = self.r_c 
            self.r_c *= (p_evap) ** (-1./3.)
            print "p_evap = %5.3f, r_c(before) = %5.3f pc, r_c(now) = %5.3f pc" % (p_evap, rc_prev/ac.pc, self.r_c/ac.pc)
        # L_c changes according the internal structure we assume
        # ----------------------------------------------------------------
        self.n_c = self.rho_c / (pc.mh * self.mu) # Assuming neutral
        p_evap = conduction.evap_pressure_fac(self.r_c, n_ps, T_ps) * p_therm_a
        self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps) * self.f_cond
        self.tau_khi = bowshock.fac_khi(self.mach) * self.tau_cc
        self.a_ram = p_ps * pi * (self.r_c) ** 2 / self.M_c / self.pressure_factor
        self.debug_flag = 1
    def dynamical_update(self, dt):
        self.eta_s = 1./conduction.fac_conduction_x(self.q_s, self.mach)
        self.tau_s = conduction.ratio_temperature(self.q_s, self.mach)
        T_ps = max(self.tau_s, 1.0)
        T_ps = T_ps * self.T_a
        n_ps = self.n_a * self.eta_s
        # The Cowie & McKee 1977 conduction rate
        self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps) * self.f_cond
        # self.tau_evap = self.M_c / mlr * self.f_cond
        # Calculate Acceleration
        p_therm_c = self.n_c * self.T_c * pc.k
        p_therm_a = self.n_a * self.T_a * pc.k
        p_ps = n_ps * T_ps * pc.k
        # self.rho_c = self.rho_c0 * p_ps / p_therm_a
        self.a_ram = p_ps * pi * self.r_c ** 2 / self.M_c #/ self.pressure_factor
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
        if(VERBOSE == True):
            print "Mc=%4.2f t(KHI,evap,cc)=%5.2f, %5.2f, %5.2f mach=%4.2f Pevap/Pth=%4.2f" % (self.M_c / params.M_c, self.tau_khi / ac.myr, self.tau_evap / ac.myr, self.tau_cc / ac.myr, self.mach, P_evap)
        # ---------------- The Expansion ----------------
        self.tclock -= dt
        self.rho_c = self.rho_c0 * (p_ps / p_therm_a / self.pressure_factor)
        if(ADJUST_RCLOUD_FOR_VAPORPRESSURE):
            self.rho_c *= P_evap
        self.r_c = (0.75 * self.M_c / self.rho_c / pi) ** (1./3.) # # SPHERICAL cloud        
        # self.r_c = sqrt(self.M_c / (self.N_c * pc.mh * self.mu) / pi) 
        self.n_c = self.rho_c / (pc.mh * self.mu)
        # if(self.time < RES_T * self.tau_cc * 3):
        #     print self.time, self.M_c, self.n_c, (self.v_h - self.vrel)/1.e5, self.L/ac.pc, self.v_max, self.r_c/ac.pc, self.a_ram, p_ps / p_therm_a

params = parameters("x300v1000")
VERBOSE = False # Output every timestep?
def run(params):
    cloud = study_case(params)
    M_0 = params.M_c
    t75, t50, t10 = -1., -1., -1.
    cloud.derived_quantities()
    cloud.cloud_shock()
    tevap = cloud.tau_evap / cloud.tau_cc
    dt = cloud.tau_cc * RES_T

    print "Initial: n_c= %5.2f; r_c= %4.1f; L_c= %4.1f" % (
        cloud.n_c, cloud.r_c / ac.pc, cloud.L / ac.kpc)
    print "    ---  Tps/T0 = %5.2f; Pps/P1 = %5.2f; v_max= %4.1f; v_i= %4.1f" % (
        cloud.tau_s, cloud.eta_s * cloud.tau_s, cloud.L/ac.kpc,
        cloud.a_ram * (cloud.r_c0 / cloud.r_c) ** 2 * cloud.tau_cc / 1.e5)

    ttot = 0.0
    # t, dv, r, R_c, L_c, n_c
    cloud.t75 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    cloud.t50 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    cloud.t25 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]        
    cloud.t10 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]    
    while ttot < END_T:
        cloud.dynamical_update(dt)
        if(cloud.M_c / M_0 < 0.75 and cloud.t75[0] == -1):
            cloud.t75 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
            print "t75/tcc= %5.1f; v= %6.3f; r_c= %4.1f; L= %4.1f; n_c= %4.1f" % (
                cloud.time / cloud.tau_cc, cloud.t75[1], cloud.t75[3], cloud.t75[4], cloud.t75[5])
        if(cloud.M_c / M_0 < 0.5 and cloud.t50[0] == -1):
            cloud.t50 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
            print "t50/tcc= %5.1f; v= %6.3f; r_c= %4.1f; L= %4.1f; n_c= %4.1f" % (
                cloud.time / cloud.tau_cc, cloud.t50[1], cloud.t50[3], cloud.t50[4], cloud.t50[5])
        if(cloud.M_c / M_0 < 0.25 and cloud.t25[0] == -1):
            cloud.t25 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
            print "t25/tcc= %5.1f; v= %6.3f; r_c= %4.1f; L= %4.1f; n_c= %4.1f" % (
                cloud.time / cloud.tau_cc, cloud.t25[1], cloud.t25[3], cloud.t25[4], cloud.t25[5])
        if(cloud.M_c / M_0 < 0.10 and cloud.t10[0] == -1):
            cloud.t10 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
            print "t10/tcc= %5.1f; v= %6.3f; r_c= %4.1f; L= %4.1f; n_c= %4.1f" % (
                cloud.time / cloud.tau_cc, cloud.t10[1], cloud.t10[3], cloud.t10[4], cloud.t10[5])
        if(cloud.M_c / M_0 < 0.1):
            print "STOP. Cloud is nearly disrupted. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
                (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
            break
        if(cloud.mach < 0.5):
            print "STOP. Flow becomes subsonic. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
                (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
            break
        ttot += RES_T
    print "End. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
        (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
    tcc = cloud.tau_cc / ac.myr
    line = "%s\t%4.2f\t" % (cloud.model, cloud.tau_cc / ac.myr)
    line += "%4.1f\t%3.1f\t" % (cloud.v_shock_approx/1.e5, 0.0)
    # line += "%4.1f\t%3.1f\t" % (cloud.v_shock_approx/1.e5, cloud.v_max0)    
    line += "%4.2f\t" % (cloud.q_s)
    line += "%4.1f\t%4.1f\t%4.1f\t" % (cloud.t75[0]/tcc, cloud.t50[0]/tcc, cloud.t25[0]/tcc)
    line += "%4.1f\t%4.1f\t%4.1f\n" % (cloud.t75[1], cloud.t50[1], cloud.t25[1])
    print line
    return cloud, line

print "Compiled."

def batch():
    fmodels = open("modellist.txt", "r")
    fout = open("models_spherical_190912.dat", "w")
    i = 0
    fout.write("#ID	Modelname	tcc	vs	vexp	q	t75	t50	t25	v75	v50	v25	Group\n")
    for modelname in fmodels:
        cld, line = run(parameters(modelname.split()[0]))
        fout.write("%d\t"%(i)+line)
        i += 1
    fout.close()
    fmodels.close()

# note 1: How do we define R_c right after the cloud shock?
# We have calculated the post-shock density rho_c and assume pi * R_c^3 * rho_c = M_c
# Now L_c is not necessarily equal to R_c. However, it is the pressure at the HEAD that determines the R_c
# The cloud will expand and L_c will increase, but the R_c will not.
# Still, it's interesting to see what a factor will be applied to R_c if we assume NOT a constant density, but a simple wave structure in the cloud:

# note 2: Expansion
# A failed tempt:
        # flux_cond = mlr_tail / (4.0 * pi * self.r_c ** 2)
        # flux_mass = self.rho_c * exp(-self.v_max) * self.v_max * self.cs_c
        # if(flux_cond <= flux_mass and self.rho_c * exp(-self.v_max) > self.rho_c0): # Only expand when evaporation is inefficient
        #     self.L = self.L + v_exp * dt
        # else:
        #     # Otherwise expansion pauses for tclock = 0.1 * t
        #     # and set the v_max to a lower value
        #     if(self.tclock <= 0):
        #         print "Lmax reached: t= %5.1f, L= %5.2f, v_max = %3.1f" % \
        #             (self.time/self.tau_cc, self.L/ac.pc, self.v_max)
        #         self.v_max = 0.9 * self.v_max
        #         if(self.v_max < 1.0):
        #             self.v_max = 1.0
        #             print "WARNING: v_exp drops below sound speed"
        #             # At this point, evaporation should destroy cloud very soon
        #         self.tclock = (1.0 - 0.9) * self.v_max * self.time
        #     # Set a timer when the new tail end reaches L
        # self.L = self.L - (self.v_h - self.vrel) * dt
        # self.L = min(self.L, self.r_c * self.v_max) # Should refine later
        # if(self.tclock <= 0):
        #     surf_dens = self.rho_c * self.L / self.v_max
        # else:
        #     surf_dens = self.rho_c * self.L / self.v_max * 0.9
        # self.r_c = sqrt(self.M_c / surf_dens / pi)

