# Some cases of interests:
import conduction
import bowshock
from astroconst import pc, ac
from numpy import exp, log, pi, sqrt
import cloud
from astropy.table import Table

tab_model = Table.read("ambient.dat", format="ascii", data_start=0)

MID = 1
q_s = 0.90
F_COND = 1.0
M_0 = 6.7e4 * ac.msolar
TEQ = 1.e4
v_h = tab_model[MID-1]["v_h"] * 1.e5
T_a = tab_model[MID-1]["T_a"]
mach = v_h / bowshock.cs_iso(T_a)
par = cloud.parameters(q_s, mach)
par.M_c = M_0
par.T_eq = TEQ
par.T_a = T_a
par.v_h = v_h
par.n_a = tab_model[MID-1]["n_a"]

MCLST = [1.e5, 1.e4, 1.e3]
TEQLST = [1.e4, 1.e3, 1.e2]
outname = "model"+str(MID)+".tab"
fout = open(outname, "w")
line = "#M_c    T_eq      t75   r75   v75   n75  Rc75   Lc75   t50   r50   v50   n50  Rc50   Lc50 v_shock\n"
fout.write(line)

for mc in MCLST:
    for Teq in TEQLST:
        par.M_c = mc * ac.msolar
        par.T_eq = Teq
        par.n_c = par.n_a * par.T_a / par.T_eq
        par.f_cond = F_COND
        cld = cloud.run(par)
        line = "%3.1e %3.1e %5.1f %5.1f %5.1f %5.2f %5.1f %6.1f %5.1f %5.1f %5.1f %5.2f %5.1f %6.1f %5.1f\n" % \
               (mc, Teq, \
                cld.t75[0], cld.t75[2], cld.t75[1], cld.t75[5], cld.t75[3], cld.t75[4], \
                cld.t50[0], cld.t50[2], cld.t50[1], cld.t50[5], cld.t50[3], cld.t50[4], \
                cld.v_shock / 1.e5)
        fout.write(line)
fout.close()
