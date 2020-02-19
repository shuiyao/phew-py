# Calculate cloud conduction in a hot ambient medium
# Reference: Cowie & McKee (1977)

# Everything is carried out in c.g.s unit system

from astroconst import pc, ac
from numpy import exp, log, pi, sqrt, linspace

print "Compiling ... bowshock.py ... DONE"

GAMMA = 5./3.
GAMMAP1 = GAMMA + 1.0
GAMMAM1 = GAMMA - 1.0

# The Rankine-Hugoniot Jump Conditions:
# --------------------------------
def rh_density_ratio(mach):
    mssq = mach * mach
    return (GAMMAP1 * mssq) / (GAMMAP1 + GAMMAM1 * (mssq - 1.0))

def rh_velocity_ratio(mach):
    return 1./rh_density_ratio(mach)

def rh_mssq(mach): # post-shock mach number
    mssq = mach * mach
    mssq2 = (2. + GAMMAM1 * mssq) / (2. * GAMMA * mssq - GAMMAM1)
    return mssq2

def rh_pressure_ratio(mach):
    mssq = mach * mach
    return (GAMMAP1 + 2.0*GAMMA*(mssq - 1.0)) / GAMMAP1

def rh_temperature_ratio(mach):
    mssq = mach * mach
    y = GAMMAP1 + 2.0 * GAMMA * (mssq - 1.0)
    y *= (GAMMAP1 + GAMMAM1 * (mssq - 1.0))
    y /= GAMMAP1 ** 2 * mssq
    return y

# Ram pressure
# --------------------------------
# Subsonic:
#   - Pram = fac * rho_a * v_rel ^ 2
#   - P2/P1 = (1.0 + 0.5 * (gamma - 1) * M^2) ^ (gamma / gamma_minus1)
# Supersonic: 
#   - P2/P1 is more complicated (See MC75). 
def postshock_pressure_fac(mach):
    mssq = mach * mach
    if(mssq >= 1): # Supersonic
        y = (0.5 * GAMMAP1) ** (GAMMAP1/GAMMAM1)
        y *= (GAMMA - 0.5 * GAMMAM1 / mssq) ** (- 1.0 / GAMMAM1)
        return y * mssq
    else:
        y = (1.0 - 0.5 * GAMMAM1 * mssq) ** (GAMMA / GAMMAM1)
        return y
# Typical Values
# Mach:   ~0.0  0.5     1.0     2.0     5.0
# Factor: ~1.0  1.22    2.05    6.35    37.17

def bernouli_pressure(mach):
    # enhancement from post-shock pressure to the stagnation point
    # mach number is the pre-shock mach number
    mssq2 = rh_mssq(mach)
    return (1. - 0.5 * GAMMAM1 * mssq2) ** (- GAMMA / GAMMAM1)

def bernouli_density(mach):
    # enhancement from post-shock pressure to the stagnation point
    # mach number is the pre-shock mach number
    mssq2 = rh_mssq(mach)
    return (1. - 0.5 * GAMMAM1 * mssq2) ** (-1.0 / GAMMAM1)
    
def ram_pressure_fac(mach):
    mssq = mach * mach
    y = postshock_pressure_fac(mach)
    return (y - 1.0) / (GAMMA * mssq)
# Typical Values
# Mach:   ~0.0  0.1     0.25    0.5     0.75    0.9     1.0     2.0     5.0     ~inf
# Factor: ~0.5  0.501   0.508   0.532   0.572   0.606   0.632   0.802   0.868   0.881

def cs_iso(T):
    if(T > 8.e3): return sqrt(GAMMA * pc.k * T / 0.62 / pc.mh)
    else: return sqrt(GAMMA * pc.k * T / 1.30 / pc.mh)    

# Cloud crushing time
def tau_cloud_crushing(r_c, n_c, n_h, vrel):
    ksi = sqrt(n_c / n_h)
    return ksi * r_c / vrel

# Disruption timescale due to Kelvin-Helmholtz Instability
# Adopted from Scannapieco & Bruggen 2015, Eq. (14), or Eq. (22)
def fac_khi(mach, fkhi = 10.0):
    # return fkhi * sqrt(1. + 4.*(GAMMA - 1.)*mach*mach)
    return fkhi * sqrt(1. + mach)

# The shock velocity
# An approximation is v_shock = R_c / tau_cc
# But a better treatment is to assume p_c = p_ps after cloud shock
# In this case: v_shock = sqrt(p_ps / p_c) * cs_c
def shock_velocity(): 
    return 0

def t_life(mach, alpha=12.5):
    msq = mach * mach
    t_shock = 2.0 * mach / (msq + 1) ** (1./2.)
    t_exp = alpha * mach / (msq + 1) ** (1./3.)
    print "Shock Time: %5.3f tcc" % (t_shock)
    print "Expansion Time: %5.3f tcc" % (t_exp)
    print "Lifetime: %5.3f tcc" % (t_shock + t_exp)
    return t_shock + t_exp

