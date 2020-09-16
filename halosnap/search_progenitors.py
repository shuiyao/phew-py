from numpy import log10, linspace, sqrt, pi
from cosmology import tcosmic
from astroconst import pc, ac
import ioformat

logz = linspace(0., 1., 401)
redz = 10. ** logz - 1.
ascales = 1. / (redz[::-1] + 1.)
snapbeg, snapend = 1, 400

HPARAM = 0.7
BOXSIZE = 12000.
unit_Time = sqrt(8.0*pi/3.0)*ac.mpc/(100.e5*HPARAM)
unit_Length = BOXSIZE * ac.kpc
unit_Velocity = unit_Length / unit_Time
unit_mass = 4.3e5 * (BOXSIZE / 25000.) ** 3
XTOL = 0.02
MTOL = 1.0
fbase = "/scratch/shuiyao/data/"
model = "l12n144-phew-movie"
gid = 1477

foutname = "progen_"+("00000"+str(gid))[-5:]
fout = open(foutname, "w")

snapi = snapend
dtime = tcosmic(ascales[snapi]) - tcosmic(ascales[snapi-1])
dtime = dtime * ac.yr / unit_Time
snapstr = ("000" + str(snapi))[-3:]
fgalname = fbase + model + "/gal_z" + snapstr + ".stat"
mgal, xgal, ygal, zgal, vxgal, vygal, vzgal = ioformat.rcol(fgalname, [4, 18, 19, 20, 15, 16, 17])
idx = gid - 1
xprev = xgal[idx] - vxgal[idx] * dtime
yprev = ygal[idx] - vygal[idx] * dtime
zprev = zgal[idx] - vzgal[idx] * dtime
mprev = mgal[idx]
fout.write("%5d %7.5f %7.5f %7.5f %6.3f\n" % (idx, xprev, yprev, zprev,\
            log10(mprev * unit_mass * 1.e10 / 0.7)))

mgalz0 = log10(mprev * unit_mass * 1.e10 / 0.7)

for snapi in range(snapbeg, snapend)[::-1]:
    dtime = tcosmic(ascales[snapi]) - tcosmic(ascales[snapi-1])
    dtime = dtime * ac.yr / unit_Time
    snapstr = ("000" + str(snapi))[-3:]
    fgalname = fbase + model + "/gal_z" + snapstr + ".stat"
    mgal, xgal, ygal, zgal, vxgal, vygal, vzgal = ioformat.rcol(fgalname, [4, 12, 13, 14, 15, 16, 17])
    rmin = 1.0
    idx = -1
    for i in range(len(xgal)):
        if(abs(xgal[i] - xprev) > XTOL): continue
        if(abs(ygal[i] - yprev) > XTOL): continue
        if(abs(zgal[i] - zprev) > XTOL): continue
        if(abs(log10(mgal[i]) - log10(mprev)) > MTOL): continue
        r = sqrt((xgal[i] - xprev) ** 2 + (ygal[i] - yprev) ** 2 + (zgal[i] - zprev) ** 2)
        if(r < rmin):
            rmin = r
            idx = i
    if(idx == -1): break
    dtime = 0.0
    xprev = xgal[idx] - vxgal[idx] * dtime
    yprev = ygal[idx] - vygal[idx] * dtime
    zprev = zgal[idx] - vzgal[idx] * dtime
    mprev = mgal[idx]
    fout.write("%5d %7.5f %7.5f %7.5f %6.3f\n" % (idx, xprev, yprev, zprev,\
                log10(mprev * unit_mass * 1.e10 / 0.7)))
    print snapi, rmin
fout.close()

def show():
    import matplotlib.pyplot as plt
    x, y, z = ioformat.rcol(foutname, [1,2,3])
    plt.plot(x, z, "b.-")
    for i in range(len(x)):
        plt.text(x[i], z[i], str(400-i), fontsize=6)
    plt.axis([-0.5, 0.5, -0.5, 0.5])
    plt.show()

print mgalz0
