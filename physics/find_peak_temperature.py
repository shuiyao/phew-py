import conduction
from astroconst import pc, ac
from numpy import sqrt, exp, log

n_0 = 1. / 300.
rho_0 = n_0 * pc.mh * 0.6
v_0 = 1.7e3 * 1.e5 # [km/s]
T_0 = 3.e6
M_0 = 6.4
fac = pc.k / (0.6 * pc.mh)

def flux_energy(rho, v, T, q):
    F = 0.5 * rho * v ** 3
    F += 2.5 * rho * v * fac * T
    F += q
    return F

def flux_kinetic_energy(rho, v):
    return 0.5 * rho * v ** 3

def q_sat(n_e, T_e):
    return 0.34 * n_e * pc.k * T_e * sqrt(pc.k * T_e / pc.me)

F_0 = flux_energy(rho_0, v_0, T_0, 0.0)
Fkin_0 = flux_kinetic_energy(rho_0, v_0)
q_0 = q_sat(n_0, T_0)
print flux_kinetic_energy(rho_0, v_0) / F_0
print q_0 / F_0

q = 0.80
rho_1 = rho_0 * conduction.ratio_density(q, M_0)
v_1 = v_0 * rho_0 / rho_1
T_1 = T_0 * conduction.ratio_temperature(q, M_0)
n_1 = rho_1 / (0.6 * pc.mh)
T_sat = T_1 * 0.35
q_sat_1 = q_sat(n_1, T_sat)
F_1 = flux_energy(rho_1, v_1, T_sat, q_sat_1)
T_peak = T_1 * 1.2
F_2 = flux_energy(rho_1, v_1, T_peak, 0.0)
print F_2 / F_0
