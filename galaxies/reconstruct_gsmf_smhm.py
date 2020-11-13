from numpy import genfromtxt, log10, exp, array, linspace
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import ioformat
import behroozi

REDSHIFT = 0.0
if(REDSHIFT == 1.0):
    moster18 = "/scratch/shuiyao/sci/REFERENCES/moster_2013_2018/moster18_z1.dat"
    hmf = "/home/shuiyao/code/Python/hmfz1.dat"
if(REDSHIFT == 0.0):
    moster18 = "/scratch/shuiyao/sci/REFERENCES/moster_2013_2018/moster18_z0.dat"
    hmf = "/home/shuiyao/code/Python/hmfz0.dat"

# mh, ms = genfromtxt(moster18, unpack=True)
mh, dm, e1, e2 = behroozi.read_smmr(REDSHIFT)
mh = array(mh)
ms = mh + array(dm)
fm18 = interp1d(mh, ms)

mh, phi = genfromtxt(hmf, unpack=True)

mhb, phib = [], []
for i in range(len(mh)):
    if(11.0 < mh[i] < 13.5):
        mhb.append(mh[i]-log10(0.7))
        phib.append(log10(phi[i]/0.7**3))
        # phib.append(log10(phi[i]))
ms = fm18(mhb)

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.plot(ms, phib, "k--")

if(REDSHIFT == 1.0):
    gsmfs_tomczak14 = [ \
        "../REFERENCES/tomczak14/t14_z0.2_z0.5.dat",\
        "../REFERENCES/tomczak14/t14c_z0.75_z1.25.dat",\
        "../REFERENCES/tomczak14/t14c_z1.5_z2.5.dat",\
    ]
    colors = ["green", "red", "blue"]
    i = 1
    ms, phi, e1, e2 = ioformat.rcol(gsmfs_tomczak14[i], [0,1,2,3], linestart=3)
    lower = -array(e2)
    upper = array(e1)
    ax.errorbar(ms, phi, yerr=[lower, upper], color=colors[i], fmt='o')

if(REDSHIFT == 0.0):
    gsmf_baldry12 = "../REFERENCES/Baldry12/gsmf_baldry12_z0.dat"
    ms, phi, err = ioformat.rcol(gsmf_baldry12, [0,2,3], linestart=3)
    phi = array(phi) / 1.e3
    err= array(err) / 1.e3
    lower = log10(phi) - log10(phi - err)
    upper = log10(phi + err) - log10(phi)
    phi = log10(phi)
    ax.errorbar(ms, phi, yerr=[lower, upper], color="black", fmt='o')

plt.show()
