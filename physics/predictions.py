# Model from Bruggen and Scannapieco 2016
from numpy import sqrt, log10, log, exp
from numpy import genfromtxt

const_A = 0.01
const_etac = 0.5
const_Tevap = 3.e6
const_cevap = 263.0
gamma_minus1 = 2./3.
gamma_plus1 = 8./3.
gamma = 5./3.
const_etashock = 0.4
const_etarocket = 0.6

chi0 = 300.0
Mach = 6.46
Nc = 1.5e20
vh = 1700.0
# mlr defined as (dm/dt)/m * tcc
def find_tevap(chi0, Mach, Nc):
    g = 3.5 * (const_etac / 0.5) * (const_A / 0.01) * (3.e6 / const_Tevap)
    g = g * (Nc / 3.e20) * Mach * sqrt(1000. / chi0)
    mssq = Mach * Mach
    fM = ((gamma_minus1) * mssq + 2.0) * (2.0 * gamma * mssq - gamma_minus1)
    fM = fM / (4.0 * gamma_plus1 * gamma_plus1 * mssq)
    fM = max(1.0, fM)
    mlr = const_A * fM * sqrt(chi0) * (1.0 - sqrt(1.0 + 4.0 * g)) / (2.0 * g)
    return 1./(-mlr)

def find_vc(t_evap, vh):
    v1 = const_etashock * vh / sqrt(chi0) * sqrt(Mach / 30. * t_evap)
    v2 = const_etarocket * const_cevap * sqrt(Mach / 30. / t_evap)
    return v1 + v2

tab = genfromtxt("models_bs16_physical.dat", names=True, dtype=('i8,S20,f8,f8,f8,f8,f8,f8'))

fout = open("models_bs16_predictions.dat", "w")
fout.write("#ID	Modelname	tevap	vc\n")
for i in range(len(tab)):
    t_evap = find_tevap(tab['chi0'][i], tab['M'][i], tab['ncRc'][i])
    vc = find_vc(t_evap, tab['vh'][i])
    line = "%d       %12s     %4.1f   %4.1f\n" % (i, tab['Modelname'][i], t_evap, vc)
    fout.write(line)
fout.close()
