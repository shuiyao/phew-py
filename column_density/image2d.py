from myinit import *
import matplotlib.pyplot as plt

ionnames = ["XX","OVI","CIV","NV","CII","NeVIII","CIII","MgII","SiIII","SiIV","HI"]
fbase = "/home/shuiyao_umass_edu/scidata/Jneil/column_density/otherBackgrounds/"
# fbase = "/home/shuiyao_umass_edu/scidata/Jneil/column_density/z1.25/"
# model = "LowCond_v1700_chi300_cond"
# model = "HC_v1000_chi300_cond"

model1 = "T0.3_v1700_chi300_cond_1e6"
model2 = "T0.3_v1700_chi300_cond_1e6T2"

fig, axs = plt.subplots(2, 2, figsize=(7,7))

def show2d(ax, model, ionid):
    NPIX = 640000
    fnamex = fbase + model + "_x.csv"
    fnamey = fbase + model + "_y.csv"
    #cols = [20+ion]
    cols = [ionid]
    ioncds = ioformat.rcol(fnamey, cols, separator=",", linestart=1)
    ioncds = array(ioncds).reshape(800, 800)
    ioncds = log10(ioncds)
    ax.pcolor(ioncds, cmap=plt.get_cmap("jet"), vmin=12., vmax=16.5)
    ax.set_xlim(300., 500.)
    ax.set_ylim(300., 500.)    

ionid = 5 # NeVIII
ionid = 10 # HI
ionname = ionnames[ionid]
show2d(axs[0][0], model1, ionid)
axs[0][0].set_title(ionname+", T = 1e4")
show2d(axs[0][1], model2, ionid)
axs[0][1].set_title(ionname+", T = 2e4")
ionid = 1 # OVI
ionname = ionnames[ionid]
show2d(axs[1][0], model1, ionid)
axs[1][0].set_title(ionname+", T = 1e4")
show2d(axs[1][1], model2, ionid)
axs[1][1].set_title(ionname+", T = 2e4")
fig.subplots_adjust(left=0.1, right=0.9, wspace=0.01, hspace=0.1)
from pylab import setp
for ax in axs.reshape(4):
    setp(ax.get_xticklabels(), visible=False)
    setp(ax.get_yticklabels(), visible=False)    
plt.savefig("/home/shuiyao_umass_edu/figures/tmp.png")
plt.show()

print ("done.")
