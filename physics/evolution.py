# Some cases of interests:
import conduction
import bowshock
from astroconst import pc, ac
from numpy import exp, log, pi, sqrt, logspace, log10
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['mathtext.default'] = "regular"
mpl.rcParams['axes.labelsize'] = "large"

class parameters():
    def __init__(self, kind=0, q=0.0, mach=0.0):
        self.ram = 1.0
        self.p_evap = 1.0
        self.t_evap = 1.0
        self.n_ps = 1.0
        self.T_ps = 1.0
        if(kind == 1):
            self.n_ps = conduction.fac_conduction_density(q, mach)
            self.T_ps = conduction.fac_conduction_temperature(q, mach)

class study_case():
    def __init__(self, M_c, T_eq, vrel):
        self.M_c = M_c
        self.vrel = vrel
        self.v_h = vrel
        self.T_c = T_eq
        self.T_a = 3.e6
        self.n_c = 1.0
        self.n_a = self.n_c / 300.0
        self.v_c = 0.0
        self.a_ram = 0.0
        self.r = 0.0
        self.time = 0.0
    def derived_quantities(self, verbose=False):
        self.rho_c = self.n_c * 0.62 * pc.mh
        self.rho_a = self.n_a * 0.62 * pc.mh        
        self.r_c = (self.M_c / (4.18879 * self.rho_c)) ** (1./3.)
        self.cs_a = bowshock.cs_iso(self.T_a)
        self.cs_c = bowshock.cs_iso(self.T_c)
        self.p_c = self.n_c * self.T_c * pc.k
        self.mach = self.vrel / self.cs_a
        self.tau_cc = bowshock.tau_cloud_crushing(self.r_c, self.n_c, self.n_a, self.vrel)
        self.N_c = self.r_c * self.n_c
        if(verbose == True):
            print "Initial Conditions"
            print "--------------------------------"
            print "R_c = ", self.r_c / ac.kpc * 1.e3, "[pc]"
            print "cs_c = ", self.cs_c / 1.e5, "[km/s]"
            print "cs_a = ", self.cs_a / 1.e5, "[km/s]"
            print "Mach = ", self.mach
            print "N_c = ", self.N_c, "[cm^-2]"
            print "tau_cc = ", self.tau_cc / 1.e6 / ac.yr, "[Myr]"
    def cloud_shock(self, conduction_flag=False, verbose=False):
        v_shock_approx = self.r_c / self.tau_cc
        fac_ram = bowshock.ram_pressure_fac(self.mach) * params.ram
        # We should use the conductive post-shock condition
        # p_ps = fac_ram * self.rho_a * self.vrel ** 2
        if(conduction_flag == True):
            T_ps = bowshock.rh_temperature_ratio(self.mach) * params.T_ps
            if(T_ps < 1): T_ps = 1
            T_ps = T_ps * self.T_a
            n_ps = self.n_a * bowshock.rh_density_ratio(self.mach) * params.n_ps
            # self.rho_c *= conduction.evap_pressure_fac(self.r_c, self.n_a, self.T_a)
        else:
            n_ps, T_ps = self.n_a, self.T_a
        p_ps = n_ps * T_ps * pc.k
        p_therm_c = self.n_c * self.T_c * pc.k
        p_therm_a = self.n_a * self.T_a * pc.k
        print "Pram Vs. P_ps = ", log10(fac_ram * self.rho_a * self.vrel ** 2), log10(p_ps)
        v_shock = sqrt(p_ps / self.p_c) * self.cs_c
        # self.rho_c = self.rho_a * (self.vrel / self.cs_c) ** 2
        if(self.mach < 1.2):
            p_ps *= bowshock.postshock_pressure_fac(self.mach)
        self.rho_c = self.rho_c * p_ps / p_therm_a
        # self.rho_c = self.rho_c * (v_shock / self.cs_c) ** 2
        self.r_c = (self.M_c / self.rho_c / pi) ** (1./3.) # cylindrical geometry
        # self.rho_c *= conduction.evap_pressure_fac(self.r_c, n_ps, T_ps)
        # print "fac =", conduction.evap_pressure_fac(self.r_c, n_ps, T_ps)
        self.n_c = self.rho_c / (pc.mh * 0.62)
        print "c_evap = ", conduction.evap_pressure_fac(self.r_c, n_ps, T_ps) * p_therm_a * \
            4. * pi * self.r_c ** 2 / conduction.mass_loss_rate(self.r_c, n_ps, T_ps) / (1.e5)
        self.r_c = (self.M_c / self.rho_c / pi) ** (1./3.) # cylindrical geometry
        v_ps = v_shock - self.cs_c ** 2 / v_shock # post-shock velocity in rest-frame
        p_cond = conduction.evap_pressure_fac(self.r_c, n_ps, T_ps) * p_therm_a
        self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps)
        self.tau_khi = bowshock.fac_khi(self.mach) * self.tau_cc        
        self.N_c = self.M_c / (0.62 * pc.mh * pi * self.r_c * self.r_c)
        self.L = self.r_c
        # self.vrel = self.vrel - v_shock
        fac_ram = bowshock.ram_pressure_fac(self.mach) * params.ram
        p_ps = fac_ram * self.rho_a * self.vrel ** 2
        # self.a_ram = (p_ps - p_therm_a) * pi * self.r_c ** 2 / self.M_c
        self.a_ram = (p_ps - self.T_a * self.n_a * pc.k) * pi * (0.1*ac.kpc) ** 2 / self.M_c        
        # self.vrel = self.vrel - self.a_ram * self.tau_cc
        if(verbose == True):
            print ""
            print "Post Cloud Shock"
            print "--------------------------------"
            print "n_c = ", self.n_c, "[cm^-3]"        
            print "R_c = ", self.r_c / ac.kpc * 1.e3, "[pc]"
            print "v_shock = ", v_shock / 1.e5, "[km/s]"
            print "v_shock_approx = ", v_shock_approx / 1.e5, "[km/s]"                
            print "Pc/k = ", p_therm_c / pc.k, "[K]"
            print "Pa/k = ", p_therm_a / pc.k, "[K]"
            print "Pps/k = ", p_ps / pc.k, "[K]"
            print "Pcond/k = ", p_cond / pc.k, "[K]"
            print "sigma0 = ", conduction.sigma_cond(self.r_c, n_ps, T_ps)
            print "t_cc = ", self.tau_cc / 1.e6 / ac.yr, "[Myr]"
            print "t_khi = ", self.tau_khi / 1.e6 / ac.yr, "[Myr]"
            print "t_evap = ", self.tau_evap / 1.e6 / ac.yr, "[Myr]"
            print "N_c = ", self.N_c, "[cm^-2]"        
    def dynamical_update(self, dt, verbose=False, conduction_flag=True, compression_flag=False, geometry="spherical"):
        if(conduction_flag == True):
            # TRICK
            # self.a_ram *= conduction.evap_pressure_fac(self.r_c, self.n_a, self.T_a)
            T_ps = bowshock.rh_temperature_ratio(self.mach) * params.T_ps
            if(T_ps < 1): T_ps = 1
            T_ps = T_ps * self.T_a
            n_ps = self.n_a * bowshock.rh_density_ratio(self.mach) * params.n_ps
        else:
            T_ps, n_ps = self.T_a, self.n_a
        self.tau_evap = conduction.tau_evap(self.M_c, self.r_c, n_ps, T_ps)
        # print self.tau_evap / ac.yr / 1.e6, n_ps, T_ps, self.r_c, self.L
        if(geometry == "cylindrical"):
            self.tau_evap = self.tau_evap / (self.L / (2.0 * self.r_c))
            # self.tau_evap = self.tau_evap
        mlr = conduction.mass_loss_rate(self.r_c, n_ps, T_ps) * 0.5
        # mlr += conduction.mass_loss_rate(self.r_c, FACDENS*self.n_a, FACTEMP*self.T_a) * self.L / (2.0 * self.r_c)
        mlr += conduction.mass_loss_rate(self.r_c, n_ps, T_ps) * self.L / (2.0 * self.r_c)
       # mlr += conduction.mass_loss_rate(self.r_c, self.n_a, self.T_a) * self.L / (2.0 * self.r_c)
        self.tau_evap = self.M_c / mlr
        # print "MLR[Total] / MLR[Head] = ", mlr / (conduction.mass_loss_rate(self.r_c, n_ps, T_ps) * 0.5)
        # Calculate Acceleration
        fac_ram = bowshock.ram_pressure_fac(self.mach) * params.ram
        p_therm_c = self.n_c * self.T_c * pc.k
        # p_therm_a = n_ps * T_ps * pc.k
        p_therm_a = self.n_a * self.T_a * pc.k
        # print "c_evap = ", conduction.evap_pressure_fac(self.r_c, self.n_a, self.T_a) * p_therm_a * \
        #     4. * pi * self.r_c ** 2 / conduction.mass_loss_rate(self.r_c, self.n_a, self.T_a) / (1.e5)
        # p_ps = fac_ram * self.rho_a * self.vrel ** 2
        p_ps = n_ps * T_ps * pc.k
        self.a_ram = (p_ps - p_therm_a) * pi * self.r_c ** 2 / self.M_c
        print "Log(P_ram) = ", log10(p_ps)
        self.tau_khi = bowshock.fac_khi(self.mach) * self.tau_cc
        # Update Cloud Properties
        self.time = self.time + dt
        self.r = self.r + self.vrel * dt
        self.vrel = self.vrel - self.a_ram * dt
        self.mach = self.vrel / self.cs_a
        self.M_c = self.M_c * exp(-dt / self.tau_evap)
        # self.M_c = self.M_c * exp(-dt / self.tau_khi)
        # TRICK
        # rho_eq = conduction.evap_pressure_fac(self.r_c, n_ps, T_ps) * (n_ps / self.n_a) * self.rho_a * T_ps / self.T_c
        # print rho_eq, self.rho_c, n_ps, T_ps, conduction.evap_pressure_fac(self.r_c, self.n_a, T_ps)
        # rho_eq = self.rho_a * self.T_a / self.T_c
        # print rho_eq, self.n_a, self.T_a, conduction.evap_pressure_fac(self.r_c, self.n_a, self.T_a)
        # rho_eq = conduction.evap_pressure_fac(self.r_c, self.n_a, self.T_a) * self.rho_a * self.T_a / self.T_c
        rho_eq = self.rho_a * self.T_a / self.T_c
        if(compression_flag == False):
            L_eq = self.M_c / (rho_eq * self.r_c * self.r_c * pi)
            # v_exp = self.cs_c * 3.0
            # v_exp = VEXP * 1.e5
            v_exp = self.cs_c * sqrt(1.0 + self.mach * self.mach) - (self.v_h - self.vrel)
            if(v_exp < 0): v_exp = 0
            self.L = self.L + v_exp * dt
            # Here using free expansion approximation: v = 2.0 * c / (gamma - 1)
            self.r_c = sqrt(self.M_c / (self.N_c * pc.mh * 0.62) / pi)
            if(self.L > L_eq): self.L = L_eq
        else:
            L_eq = self.M_c / ((self.T_a * self.rho_a / self.T_c) * self.r_c * self.r_c * pi)
            # Lower pressure (no vapor pressure), therefore longer L_eq
            self.L = self.L + self.cs_c * dt * 3.0
            if(self.L > L_eq): self.L = L_eq
            self.rho_c = self.M_c / (self.L * self.r_c * self.r_c * pi)
            # Allow cloud size to shrink
            # However, problem is it induces a positive feedback:
            # Smaller Rc leads to higher vapor pressure and higher rho_eq, the cloud will collapse!
            if(self.rho_c < rho_eq):
                self.r_c = sqrt(self.M_c / (rho_eq * self.L * pi))

        self.rho_c = self.M_c / (self.L * self.r_c * self.r_c * pi)
        self.n_c = self.rho_c / (pc.mh * 0.62)
        if(verbose == True):
            print "%5.3f(%5.3f) %6.2f %6.1f %5.3f %6.2f(%5.3f) %6.2f %5.3f %7.2f" % \
                (self.time/ac.yr/1.e6, self.time/self.tau_cc, self.r / ac.kpc, self.vrel / 1.e5, self.M_c / M_0, \
                 self.L * 1.e3 / ac.kpc, self.L / L_eq, self.r_c * 1.e3 / ac.kpc, \
                 self.n_c, self.tau_evap / ac.yr / 1.e6)

RES_T = 0.1 # time resolution, dt = tau_cc * RES_T
END_T = 200. # End time, tend = tau_cc * END_T
M_S = 1.5106        
# params = conduction.parameters(M_S) # Doesn't matter to the calculation ...
# Above are default values for conduction calculate.

#params = parameters(kind="default")

M_0 = 6.7e4 * ac.msolar

def cloud_shock(cond = True):
    v_h = 1000. # [km/s]
    T_eq = 1.e4
    m4v1000 = study_case(M_0, T_eq, v_h * 1.e5)
    m4v1000.T_a = 3.e6
    m4v1000.n_a = m4v1000.n_c / 300.
    t50, t25 = -1., -1.
    m4v1000.derived_quantities()
    m4v1000.cloud_shock(conduction_flag=cond, verbose=True)
    return m4v1000

VEXP = 35.
vhs = [1000., 1700., 3000., 1700., 3000., 3000.]
chis = [300., 300., 300., 1000., 3000., 300.]
machs = [3.8, 6.4, 11.4, 3.6, 11.4, 11.4]
FACDENS, FACTEMP = 2.75, 1.26 # x300v1000
FACDENS, FACTEMP = 3.0, 2.5 # x300v1700
FACDENS, FACTEMP = 1.0, 1.0 # x300v1700
MIDX = 1
qs = conduction.solve_for_q(machs[MIDX])
params = parameters(kind=1, q=qs, mach=machs[MIDX])
def special_case(mid = MIDX, cond = True, plotting = False):
    # params = parameters(kind=1, q=qs, mach=machs[mid])
    v_h = vhs[mid] # [km/s]
    T_eq = 1.e4
    m4v1000 = study_case(M_0, T_eq, v_h * 1.e5)
    if(mid == 5):
        m4v1000.n_c *= 10.0
    m4v1000.T_a = 3.e6 * chis[mid] / 300.
    m4v1000.n_a = m4v1000.n_c / chis[mid]
    t75, t50, t25 = -1., -1., -1.
    m4v1000.derived_quantities()
    m4v1000.cloud_shock(conduction_flag=cond, verbose=True)
    print ""
    print "Evolution: "
    print "--------------------------------"
    print "time[Myr] r[kpc]  Vrel[km/s] Mc/M0   L[pc] Rc[pc] n_c[cm^-3] t_evap [Myr]"
    print "---"
    tevap = m4v1000.tau_evap / m4v1000.tau_cc
    dt = m4v1000.tau_cc * RES_T

    ttot = 0.0
    tx, r, mc, v = [], [], [], []
    rc25, L25, dv25, rho25 = -1.0, -1.0, -1.0, -1.0
    rc50, L50, dv50, rho50 = -1.0, -1.0, -1.0, -1.0
    rc75, L75, dv75, rho75 = -1.0, -1.0, -1.0, -1.0    
    while ttot < END_T:
        m4v1000.dynamical_update(dt, conduction_flag=cond, compression_flag=False, geometry="cylindrical", verbose=True)
        r.append(m4v1000.r / ac.kpc)
        v.append(m4v1000.vrel / v_h / 1.e5)
        mc.append(m4v1000.M_c / M_0)
        tx.append(ttot)
        if(m4v1000.M_c / M_0 < 0.75 and t75 == -1):
            t75 = ttot
            rc75 = m4v1000.r_c / ac.kpc * 1.e3
            L75 = m4v1000.L / ac.kpc * 1.e3
            dv75 = v_h - m4v1000.vrel / 1.e5
            rho75 = log10(m4v1000.rho_c)
        if(m4v1000.M_c / M_0 < 0.5 and t50 == -1):
            t50 = ttot
            rc50 = m4v1000.r_c / ac.kpc * 1.e3
            L50 = m4v1000.L / ac.kpc * 1.e3
            dv50 = v_h - m4v1000.vrel / 1.e5
            rho50 = log10(m4v1000.rho_c)
        if(m4v1000.M_c / M_0 < 0.25 and t25 == -1):
            t25 = ttot
            rc25 = m4v1000.r_c / ac.kpc * 1.e3
            L25 = m4v1000.L / ac.kpc * 1.e3
            dv25 = v_h - m4v1000.vrel / 1.e5
            rho25 = log10(m4v1000.rho_c)            
        if(m4v1000.M_c / M_0 < 0.1):
            print "STOP. Cloud is nearly disrupted."
            break
        # if(m4v1000.mach < 1):
        #     print "STOP. Flow becomes subsonic."
        #     break
        ttot += RES_T
    print "q_s = %5.3f" % (qs)
    print "tau_evap = %5.3f tau_cc" % (tevap)
    print "t75 = %5.3f tau_cc; (Rc = %5.3f pc, L = %5.3f pc, dv = %5.3f km/s, rhoc = %5.3f)" % (t75, rc75, L75, dv75, rho75)    
    print "t50 = %5.3f tau_cc; (Rc = %5.3f pc, L = %5.3f pc, dv = %5.3f km/s, rhoc = %5.3f)" % (t50, rc50, L50, dv50, rho50)
    print "t25 = %5.3f tau_cc; (Rc = %5.3f pc, L = %5.3f pc, dv = %5.3f km/s, rhoc = %5.3f)" % (t25, rc25, L25, dv25, rho25)
    mach = m4v1000.v_h / m4v1000.cs_a
    vs = m4v1000.cs_c * sqrt(mach ** 2 + 1.0) / 1.e5
    print "Shock Velocity: %5.3f [km/s]" % (vs)
    if(plotting == True):
        plt.plot(tx, v, "r-")
        plt.plot(tx, mc, "b-")
        plt.show()
    return m4v1000

def x300v1000():
    v_h = 1000. # [km/s]
    cmap1 = plt.get_cmap("cool")
    cmap2 = plt.get_cmap("cool")    
    Mci = [1.e2, 1.e3, 1.e4, 6.7e4]

    ls = [":", "--", "-"]
    lines = []
    i = 0
    for T_eq in [500., 2000., 10000.]:
        clridx = 0
        for M_0 in Mci:
            clridx += 60
            M_0 = M_0 * ac.msolar
            x300v1000 = study_case(M_0, T_eq, v_h * 1.e5)
            # x300v1000.T_a = 1.e7
            # x300v1000.n_a = x300v1000.n_c / 1000.
            x300v1000.derived_quantities(verbose=False)
            x300v1000.cloud_shock(verbose=False)
            print "Parameters: Teq = %5.0f; M0 = %4.0f; tau_cc = %5.3f; tau_evap = %5.3f" % \
            (T_eq, M_0/ac.msolar, x300v1000.tau_cc/1.e6/ac.yr, x300v1000.tau_evap/1.e6/ac.yr)
            dt = x300v1000.tau_cc * RES_T

            ttot = 0.0
            r, mc, v = [], [], []
            while ttot < END_T:
                x300v1000.dynamical_update(dt, geometry="cylindrical", verbose=False)
                r.append(x300v1000.r / ac.kpc)
                v.append(x300v1000.vrel / 1.e8)
                mc.append(x300v1000.M_c / M_0)
                if(x300v1000.M_c / M_0 < 0.1):
                    print "STOP. Cloud is nearly disrupted. t = %f tau_cc; v = %f; L/Rc = %f" % (ttot, x300v1000.vrel/1.e5, x300v1000.L / x300v1000.r_c)
                    break
                if(x300v1000.mach < 1):
                    print "STOP. Flow becomes subsonic., t = %f tau_cc; v = %f" % (ttot, x300v1000.vrel/1.e5)
                    break
                ttot += RES_T
            p, = plt.plot(r, mc, color=cmap1(clridx), linestyle=ls[i], linewidth=2)
            # plt.plot(r, v, color=cmap2(clridx), linestyle=ls[i], linewidth=2)            
            lines.append(p)
            if(len(lines) == 12):
                plt.plot(r, mc, "k-", linewidth=2)
        i = i + 1
    plt.title(r'$\chi_0 = 300; T_a = 3\times10^6 [K]; n_a = 1.6\times10^{-3} [cm^-2]$')
    plt.xlabel(r'$r [kpc]$')
    plt.ylabel(r'$\frac{M_c(t)}{M_c(0)}$', rotation=0)
    l1 = plt.legend([lines[8], lines[9], lines[10], lines[11]], [r'$10^2M_\odot$', r'$10^3M_\odot$', r'$10^4M_\odot$', r'$6.7\times10^4M_\odot$'], fontsize=12)
    plt.gca().add_artist(l1)
    plt.legend([lines[0], lines[4], lines[8]], [r'$100 K$', r'$1000 K$', r'$10000 K$'], fontsize=12, loc=4)    
    plt.axis([0., 150., 0., 1.])
    plt.show()

def x300v1000scatter():
    fig = plt.figure(1, figsize=(6, 6))
    v_h = 1700. # [km/s]
    cmap1 = plt.get_cmap("autumn")
    cmap2 = plt.get_cmap("cool")    
    Mci = logspace(2., 5., 10)
    END_T = 250

    ls = [":", "--", "-"]
    mks = ["o", "+", "^"]
    lines = []
    i = 0
    for T_eq in [500., 2000., 10000.]:
        clridx = 0
        for M_0 in Mci:
            clridx += 25
            print "Parameters: Teq = %5.0f; M0 = %4.0f" % (T_eq, M_0)
            M_0 = M_0 * ac.msolar
            x300v1000 = study_case(M_0, T_eq, v_h * 1.e5)
            x300v1000.derived_quantities(verbose=False)
            x300v1000.cloud_shock(verbose=False)
            dt = x300v1000.tau_cc * RES_T

            ttot = 0.0
            while ttot < END_T:
                x300v1000.dynamical_update(dt, geometry="cylindrical", verbose=False)
                r = x300v1000.r / ac.kpc
                v = x300v1000.vrel / v_h / 1.e5
                rl = x300v1000.r_c / x300v1000.L
                if(x300v1000.M_c / M_0 < 0.5): break
                if(x300v1000.mach < 1): break
                ttot += RES_T
            plt.plot(r, v, color=cmap1(clridx), marker=mks[i], markersize=6)
            plt.plot(r, rl, color=cmap2(clridx), marker=mks[i], markersize=6)            
        i = i + 1
    plt.title(r'$\chi_0 = 500; T_a = 3\times10^6 [K]; n_c = 1.6\times10^{-3} [cm^-2]$')
    plt.xlabel(r'$r [kpc]$')
    plt.ylabel(r'$\frac{y(t)}{x(t)}$', rotation=0)
    symbols = []
    for i in range(3):
        sym, = plt.plot(0., 0., markersize=6, marker=mks[i], color="black")
        symbols.append(sym)
    plt.legend(symbols, [r'$100 K$', r'$1000 K$', r'$10000 K$'], fontsize=12, loc=4)
    plt.text(5., 0.55, r'$\frac{v_{50}}{v_{i}}$', fontsize=20, color="black")
    plt.text(5., 0.45, r'$\frac{L_{50}}{r_{c;50}}$', fontsize=20, color="black")    
    plt.axis([0., 150., 0., 1.])
    norm1 = mpl.colors.Normalize(vmin=2, vmax=5)
    axcbar = fig.add_axes([0.25,0.55,0.55,0.015])
    cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=cmap1, norm=norm1, orientation="horizontal")
    cdcbar.set_ticks([2, 3, 4, 5])
    cdcbar.set_ticklabels(["2","3","4","5"])
    axcbar = fig.add_axes([0.25,0.45,0.55,0.015])
    cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=cmap2, norm=norm1, orientation="horizontal")
    cdcbar.set_ticks([2, 3, 4, 5])
    cdcbar.set_ticklabels(["2","3","4","5"])
    cdcbar.set_label(r'$\log(M_{c;i}/M_\odot)$')
    plt.show()

print "Compiled."
