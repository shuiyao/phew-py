# Model from Bruggen and Scannapieco 2016
from numpy import sqrt, log10, log, exp

const_A = 0.01
const_etac = 0.5
const_Tevap = 3.e6
gamma_minus1 = 2./3.
gamma_plus1 = 8./3.
gamma = 5./3.

chi0 = 300.0
Mach = 6.46
Nc = 1.5e20
# mlr defined as (dm/dt)/m * tcc
g = 3.5 * (const_etac / 0.5) * (const_A / 0.01) * (3.e6 / const_Tevap)
g = g * (Nc / 3.e20) * Mach * sqrt(1000. / chi0)
mssq = Mach * Mach
fM = ((gamma_minus1) * mssq + 2.0) * (2.0 * gamma * mssq - gamma_minus1)
fM = fM / (4.0 * gamma_plus1 * gamma_plus1 * mssq)
fM = max(1.0, fM)
mlr = const_A * fM * sqrt(chi0) * (1.0 - sqrt(1.0 + 4.0 * g)) / (2.0 * g)
