# Version 3.1: 191110
# Modify the vmax calculation
# Version 3.0: 191104
# Modify the conduction zone calculation
# Version 2.0:
# 1. Approximate the internal cloud as a similar wave parameterized by a polytropic index n
#    - 1.1. The index n is constrained from the expansion speed, which is 2/(n-1)*cs_c. n ~ 
# Version 1.0:
# The first version of cynlindrical cloud.
import conduction
import bowshock
import simplewave
from astroconst import pc, ac
from numpy import exp, log, pi, sqrt, logspace, log10
import matplotlib.pyplot as plt
import os.path
import sys
import pdb

import config_mpl

# Usage:
# params = parameters("x300v1000")
# run(params, output=False):

RES_T = 0.01 # time resolution, dt = tau_cc * RES_T
END_T = 200. # End time, tend = tau_cc * END_T
#M_S = 1.5106 # PHI_S = 1.1
M_S = 1.4324 # PHI_S = 1.0
FACDENS, FACTEMP = 1.0, 1.0

USE_MINIMUM_EXPANSION_SPEED = False # v_shock - 3.0 c_s ?
ADJUST_RHOC_WITH_PEVAP = False
KHI_MODIFIED_TIMESCALE = False # Modified tau_kh by exp(Rc/lambda_kh)

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
        self.n_ps = 0.0
        self.T_ps = 0.0

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
        self.lkh = 0.0
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
        # WARNING: Ignore the pressure difference between stagnation point and the shock front
        # if(self.mach < 1.2):
        #     p_ps *= bowshock.postshock_pressure_fac(self.mach)
        self.rho_c = self.rho_c * p_ps / p_therm_a # Isothermal shock condition
        if(USE_MINIMUM_EXPANSION_SPEED):
            self.v_max = self.v_shock_approx/self.cs_c - 3.0 # upper limit
            # self.v_max = self.v_shock_approx # upper limit
            self.v_max = min(self.v_max, log(self.rho_c/self.rho_c0))
        else:
            self.v_max = log(self.rho_c/self.rho_c0)
        print "v_shock_approx = %5.3f; v_shock = %5.3f; c*log(rho_c/rho_c0) = %5.3f" % (self.v_shock_approx / 1.e5, self.v_shock / 1.e5, log(self.rho_c/self.rho_c0)*self.cs_c / 1.e5)
        # self.v_max = log(self.rho_c/self.rho_c0)
        self.v_max0 = self.v_max
        # ----------------------------------------------------------------
        # note 1: How do we define R_c after cloud shock?
        self.r_c = (self.M_c / self.rho_c / pi / 2.) ** (1./3.) # cylindrical geometry
        p_evap_ratio = conduction.evap_pressure_fac(self.r_c, n_ps, T_ps)
        print " ---> P_evap initial: ", p_evap_ratio
        if(ADJUST_RHOC_WITH_PEVAP):
            if(p_evap_ratio > 1):
                self.rho_c *= p_evap_ratio
                self.r_c *= (p_evap_ratio) ** (-1./3.)
        # L_c changes according the internal structure we assume
        # 180712: We no longer accept the case with constant density!
        self.L = 2. * self.r_c
        self.tau_sc = 2.0 * self.r_c / self.cs_c
        if(self.n_poly == 1.0):
            # self.L = self.r_c * self.v_max # isothermal, 5.0 by default, constrained from x300v1700
            self.L = self.r_c * 4.0
            self.v_max = self.v_max
        else:
            self.L = self.r_c * (self.n_poly + 1.0) / (self.n_poly - 1.0)
        # ----------------------------------------------------------------
        self.n_c = self.rho_c / (pc.mh * self.mu) # Assuming neutral
        self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps) / self.f_cond
        self.tau_khi = bowshock.fac_khi(self.mach) * self.tau_cc
        self.N_c = self.M_c / (self.mu * pc.mh * pi * self.r_c * self.r_c)
        # self.vrel = self.vrel - v_shock
        # fac_ram = bowshock.ram_pressure_fac(self.mach) * params.ram
        # p_ps = fac_ram * self.rho_a * self.vrel ** 2
        self.a_ram = p_ps * pi * (self.r_c) ** 2 / self.M_c / self.pressure_factor
        self.debug_flag = 1
    def dynamical_update(self, dt):
        self.eta_s = 1./conduction.fac_conduction_x(self.q_s, self.mach)
        self.tau_s = conduction.ratio_temperature(self.q_s, self.mach)
        T_ps = max(self.tau_s, 1.0)
        T_ps = T_ps * self.T_a
        n_ps = self.n_a * self.eta_s
        # self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps) / self.f_cond
        # self.tau_evap = self.tau_evap / (self.L / (2.0 * self.r_c))

        # mlr = conduction.mass_loss_rate(self.r_c, n_ps, T_ps) * 0.5 # Head
        # mlr_tail = conduction.mass_loss_rate(self.r_c, self.n_a, self.T_a) 
        # mlr += mlr_tail * self.L / (2.0 * self.r_c) # Tail
        # mlr = conduction.mass_loss_rate(self.r_c, n_ps, T_ps) / self.v_max * self.L / (2.0 * self.r_c)
        # mlr = conduction.mass_loss_rate(self.r_c, n_ps, T_ps) / 3.5 * self.L / (2.0 * self.r_c)
        # mlr = 2.3e-34 * T_ps ** 2.5 * self.L / 3.5 * 2. * pi * ac.pc / 3.5

        mlr_tail = 4.46e-15 * self.T_a ** 2.5 / (self.r_c)
        self.v_max = -log(mlr_tail * self.time / self.rho_c / self.r_c)
        # SH191104: Huge bug identified here... Used to be > 0
        # SH191106: A more robust way of doing evaporation
        if(conduction.evaluate_saturation_point(self.T_c, n_ps, T_ps, self.r_c) < 0):
            # The classical conduction all the way down. No turning point.
            T_star = 1.e4
            # mlr = 4.46e-15 * (T_ps ** 2.5 - T_star ** 2.5) * self.L / 3.5 # Equation (1)
        else:
            T_star = conduction.solve_for_saturation_point(n_ps, T_ps, self.r_c, T_c=self.T_c, verbose=False)
            # mlr = 1.42e-25 * n_ps * T_ps * self.r_c * T_star ** 1.03 * self.L / 3.5 Equation (2)
            # Should be equivalent to Equation (1)
        mlr = 4.46e-15 * (T_ps ** 2.5 - T_star ** 2.5) * self.L / 3.5 # Equation (1)            

        # print n_ps, T_ps/1.e6, T_star/1.e6
            # mlr = 1.56e-15 * T_ps ** 2.5 * self.L
        self.tau_evap = self.M_c / mlr / self.f_cond
        # pdb.set_trace()
        # Calculate Acceleration
        p_therm_c = self.n_c * self.T_c * pc.k
        p_therm_a = self.n_a * self.T_a * pc.k
        p_ps = n_ps * T_ps * pc.k
        # Here we assume rho_c do not change with time???
        # if(ADJUST_RHOC_WITH_PEVAP):
        #     p_evap_ratio = conduction.evap_pressure_fac(self.r_c, n_ps, T_ps)
        #     if(p_evap_ratio > 1):
        #         self.rho_c *= p_evap_ratio
        #         self.r_c *= (p_evap_ratio) ** (-1./3.)
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
        if(KHI_MODIFIED_TIMESCALE):
            lkh = conduction.lambda_kh(1.0, self.n_c/n_ps, T_ps, n_ps, fcond=self.f_cond)
            # print "LKH, Rc = ", lkh/ac.pc, self.r_c/ac.pc
            tau_disrupt = 1./(1./self.tau_evap + exp(-lkh/self.r_c)/self.tau_khi)
            self.lkh = lkh
            self.n_ps, self.T_ps = n_ps, T_ps
        else:
            tau_disrupt = self.tau_evap
        self.M_c = self.M_c * exp(-dt / tau_disrupt)
        if(VERBOSE == True):
            print "Mc=%4.2f t(KHI,evap,cc)=%5.2f, %5.2f, %5.2f mach=%4.2f Pevap/Pth=%4.2f" % (self.M_c / params.M_c, self.tau_khi / ac.myr, self.tau_evap / ac.myr, self.tau_cc / ac.myr, self.mach, P_evap)
        # ---------------- The Expansion ----------------
        self.tclock -= dt
        self.v_max = min(self.v_max, log(self.rho_c/self.rho_c0))
        # print self.v_max, log(self.rho_c/self.rho_c0)
        # self.v_max = log(self.rho_c/self.rho_c0)
        if(self.n_poly == 1.0):
            v_exp = self.cs_c * self.v_max
        elif(self.n_poly > 1.0):
            v_exp = 2.0 * self.cs_c / (self.n_poly - 1.0)
        # flux_cond = mlr_tail
        # flux_mass = self.rho_c * exp(-self.v_max) * self.v_max * self.cs_c
        # if(flux_cond <= flux_mass and self.rho_c * exp(-self.v_max) > self.rho_c0 and self.L < 250.*ac.pc): # Only expand when evaporation is inefficient
        if(self.v_max > 0.0): # Only expand when evaporation is inefficient     
            self.L_prev = self.L
            self.L = self.L + v_exp * dt
            # if(self.time / self.tau_cc >= 15. and self.debug_flag == 1):
            #     print flux_cond / flux_mass
            #     self.debug_flag = 0
        else:
            if(self.L_prev < self.L):
                print "Lmax reached: t= %5.1f, L= %5.2f, v_max = %3.1f\n" % \
                    (self.time/self.tau_cc, self.L/ac.pc, self.v_max)
                self.L_prev = self.L
            # Set a timer when the new tail end reaches L
        # self.L = self.L - (self.v_h - self.vrel) * dt
        # ----------------------------------------------------------------
        # note 3: How does the internal structure evolve with time?
        # ----------------------------------------------------------------
        # OPTION 1: Allow rho_c (head density) to change
        # if(self.rho_c > self.rho_c0 * (p_ps / p_therm_a / self.pressure_factor)): # Density is still large
        #     self.rho_c = self.rho_c * self.L_prev / self.L
        #     self.r_c = sqrt(self.M_c / (self.rho_c * pi * self.L / 4.0))            
        # else: # Density reaches the 1/3 value
        # OPTION 2: Make rho_c strictly as constrained by the post-shock gas. Let R_c change
        self.rho_c = self.rho_c0 * (p_ps / p_therm_a / self.pressure_factor)
        # self.r_c = sqrt(self.M_c / (self.rho_c * pi * self.L / self.v_max)) # If we let the Pc = Pa
        self.r_c = sqrt(self.M_c / (self.N_c * pc.mh * self.mu) / pi) 
        self.n_c = self.rho_c / (pc.mh * self.mu)
        # if(self.time < RES_T * self.tau_cc * 3):
        #     print self.time, self.M_c, self.n_c, (self.v_h - self.vrel)/1.e5, self.L/ac.pc, self.v_max, self.r_c/ac.pc, self.a_ram, p_ps / p_therm_a

# params = parameters("x300v1700c05")
params = parameters("x1000v480")
VERBOSE = False # Output every timestep?
def run(params, q_s=None, output=False):
    cloud = study_case(params)
    M_0 = params.M_c
    t75, t50, t10 = -1., -1., -1.
    cloud.derived_quantities()
    cloud.cloud_shock()
    tevap = cloud.tau_evap / cloud.tau_cc
    dt = cloud.tau_cc * RES_T

    print "Initial: n_c= %5.2f; r_c= %4.1f; L_c= %4.1f" % (
        cloud.n_c, cloud.r_c / ac.pc, cloud.L / ac.pc)
    print "    ---  Pps/P1 = %5.2f; v_max= %4.1f; v_i= %4.1f" % (
        cloud.eta_s * cloud.tau_s, cloud.v_max,
        cloud.a_ram * (cloud.r_c0 / cloud.r_c) ** 2 * cloud.tau_cc / 1.e5)
    if(output == True):
        fout = open("cloud.dat", "w")
        fout.write("#t75/tcc Mc[g] dv[km/s] r_c[pc] L[pc] n_c[cm^-3]\n")

    if(q_s != None):
        cloud.q_s = q_s

    ttot = 0.0
    # t, dv, r, R_c, L_c, n_c
    cloud.t75 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    cloud.t50 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    cloud.t25 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]        
    cloud.t10 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]    
    while ttot < END_T:
        cloud.dynamical_update(dt)
        if(output == True):
            line = "%5.1f %7.5e %6.3f %4.1f %4.1f %4.1f\n" % (\
                cloud.time / cloud.tau_cc, cloud.M_c, \
                (cloud.v_h - cloud.vrel) / 1.e5, \
                cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                cloud.n_c)
            fout.write(line)
        if(cloud.M_c / M_0 < 0.75 and cloud.t75[0] == -1):
            cloud.t75 = [cloud.time / ac.myr, \
                         (cloud.v_h - cloud.vrel) / 1.e5, cloud.r / ac.kpc, \
                         cloud.r_c / ac.kpc * 1.e3, cloud.L / ac.kpc * 1.e3, \
                         cloud.n_c]
            print "t75/tcc= %5.1f; v= %6.3f; r_c= %4.1f; L= %4.1f; n_c= %4.1f" % (
                cloud.time / cloud.tau_cc, cloud.t75[1], cloud.t75[3], cloud.t75[4], cloud.t75[5])
            print "   - tau[evap, khi] = %4.1f %4.1f [tcc]; Lkh = %4.1f" % (
                cloud.tau_evap / cloud.tau_cc, cloud.tau_khi / cloud.tau_cc, cloud.lkh / ac.pc)
            # print "   - chi = %4.1f, n_ps = %4.2f, T_ps = %4.1e" % (
            #     cloud.n_c / cloud.n_ps, cloud.n_ps, cloud.T_ps)
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
        if(cloud.mach < 1):
            print "STOP. Flow becomes subsonic. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
                (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
            break
        ttot += RES_T
    print "End. t = %5.1f Myr; r = %5.1f kpc; v = %5.1f km/s." % \
        (cloud.time / ac.myr, cloud.r / ac.kpc, cloud.vrel / 1.e5)
    tcc = cloud.tau_cc / ac.myr
    line = "0\t%s\t%4.2f\t" % (cloud.model, cloud.tau_cc / ac.myr)
    line += "%4.1f\t%3.1f\t" % (cloud.v_shock_approx/1.e5, cloud.v_max0)
    line += "%4.2f\t" % (cloud.q_s)
    line += "%4.1f\t%4.1f\t%4.1f\t" % (cloud.t75[0]/tcc, cloud.t50[0]/tcc, cloud.t25[0]/tcc)
    line += "%4.1f\t%4.1f\t%4.1f\t" % (cloud.t75[1], cloud.t50[1], cloud.t25[1])
    print line
    return cloud

def batch(foutname, q_s=None):
    batch_models = ["x300v1000", "x300v1700", "x1000v1700", "x3000v3000", "x300v3000", \
                    "x300v3000b", "x1000v480", "x3000v860", "x300v1700c5", "x300v1700c20"]
    fout = open(foutname, "w")
    fout.write("#ID	Modelname	tcc	vs	vexp	q	t75	t50	t25	v75	v50	v25	Lmax\n")
    for modelid in range(len(batch_models)):
        params = parameters(batch_models[modelid])
        cloud = study_case(params)
        if(q_s != None):
            cloud.q_s = q_s
        M_0 = params.M_c
        t75, t50, t10 = -1., -1., -1.
        cloud.derived_quantities()
        cloud.cloud_shock()
        tevap = cloud.tau_evap / cloud.tau_cc
        dt = cloud.tau_cc * RES_T
        cloud.t75 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
        cloud.t50 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
        cloud.t25 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]        
        cloud.t10 = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]    
        print "--------> ", modelid, batch_models[modelid], cloud.f_cond

        ttot = 0.0
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
            if(cloud.M_c / M_0 < 0.25 and cloud.t25[0] == -1):
                cloud.t25 = [cloud.time / ac.myr, \
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
        tcc = cloud.tau_cc / ac.myr            
        line = "%d\t%s\t%4.2f\t" % (modelid, cloud.model, cloud.tau_cc / ac.myr)
        line += "%4.1f\t%3.1f\t" % (cloud.v_shock_approx/1.e5, cloud.v_max0)
        line += "%4.2f\t" % (cloud.q_s)
        line += "%4.1f\t%4.1f\t%4.1f\t" % (cloud.t75[0]/tcc, cloud.t50[0]/tcc, cloud.t25[0]/tcc)
        line += "%4.1f\t%4.1f\t%4.1f\t%4.0f\n" % (cloud.t75[1], cloud.t50[1], cloud.t25[1], cloud.L/ac.pc)
        fout.write(line)
    fout.close()

print "Compiled."

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

