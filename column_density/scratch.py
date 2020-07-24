import ioformat
from mymod import *
import matplotlib.pyplot as plt

# x is perpendicular to the wind

fbase = "/scratch/shuiyao/Jneil/column_density/"
# fnamex = fbase + "HC_v1000_chi300_cond_x.csv"
# fnamey = fbase + "HC_v1000_chi300_cond_x.csv"
#fnamey = fbase + "T0.3_v1700_chi300_cond_x.csv"
fnamex = fbase + "LowCond_v1700_chi300_cond_x.csv"
fnamey = fbase + "LowCond_v1700_chi300_cond_x.csv"
#fnamey = fbase + "LowCond_v1700_chi300_cond.csv"

clrs_t = ["magenta", "red", "orange", "yellow"]

snapt = "t50"

if(snapt == "t90"):
    cols = [1, 2, 10]
if(snapt == "t50"):
    cols = [21, 22, 30]    
if(snapt == "t25"):
    cols = [31, 32, 40]

# f n= open(fname, "r")
# header = f.readline()
# f.close()

cols1 = [21, 22, 30]
cols2 = [31, 32, 40]

def convert_list(lst):
    lst = array(lst)
    lst = lst[lst > 0]
    lst = log10(lst)
    return lst

# 1   2   3  4   5      6    7    8     9    10
# OVI CIV NV CII NeVIII CIII MgII SiIII SiIV HI

OVIx, CIVx, HIx = ioformat.rcol(fnamex, cols1, separator=",", linestart=1)

OVIx = convert_list(OVIx)
CIVx = convert_list(CIVx)
HIx = convert_list(HIx)

OVIy, CIVy, HIy = ioformat.rcol(fnamey, cols2, separator=",", linestart=1)

OVIy = convert_list(OVIy)
CIVy = convert_list(CIVy)
HIy = convert_list(HIy)

# step = 100
# plt.plot(CIV0[::step], OVI0[::step], "b.", alpha=0.2, markersize=4)
# plt.plot(HI2[::step], OVI2[::step], "r.", alpha=0.2, markersize=4)
# plt.plot(CIV2[::step], OVI2[::step], "r.", alpha=0.2, markersize=4)
nbins = 100
plt.hist(HI0x, bins=nbins, color="cyan", histtype="step")
plt.hist(CIV0x, bins=nbins, color="lime", histtype="step")
plt.hist(OVI0x, bins=nbins, color="orange", histtype="step")
plt.hist(HI0y, bins=nbins, color="blue", histtype="step")
plt.hist(CIV0y, bins=nbins, color="green", histtype="step")
plt.hist(OVI0y, bins=nbins, color="red", histtype="step")
# plt.hist(HI3, bins=nbins, color="cyan", histtype="step")
# plt.hist(CIV3, bins=nbins, color="lime", histtype="step")
# plt.hist(OVI3, bins=nbins, color="orange", histtype="step")
# plt.hist(HI3y, bins=nbins, color="blue", histtype="step")
# plt.hist(CIV3y, bins=nbins, color="green", histtype="step")
# plt.hist(OVI3y, bins=nbins, color="red", histtype="step")
# plt.legend(["HI90","CIV90","OVI90","HI25","CIV25","OVI25"])
# plt.legend(["HI25x","CIV25x","OVI25x","HI25y","CIV25y","OVI25y"])
plt.legend(["HI90x","CIV90x","OVI90x","HI90y","CIV90y","OVI90y"])
plt.xlabel("log(N)")
plt.ylabel("Count")
plt.yscale("log")
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

print "done"