import ioformat
from mymod import *
import matplotlib.pyplot as plt
from pltastro import frame, draw
from numpy import histogram, cumsum
import config_mpl
import os
import pdb

# x is perpendicular to the wind

PROJ = "perpendicular"
PROJ = "parallel"

fbase = "/scratch/shuiyao/Jneil/column_density/"
uvbkg = "otherBackgrounds/"
model = "T0.3_v1700_chi300_cond"
# bkgstr = ["0-01", "0-1", "10", "100"]
# bkgstr = ["1e5", "1e6", "1e7"]
bkgstr = ["1e5T2", "1e6T2", "1e7T2"]

# I. Generate the CPDF Files
for bki in [0,1,2]:
    modelstr = model + "_" + bkgstr[bki]
    fnamex = fbase + uvbkg + modelstr + "_x.csv"
    fnamey = fbase + uvbkg + modelstr + "_y.csv"
    generate_pdf_tables()

# fnamex = fbase + model + "_x.csv"
# fnamey = fbase + model + ".csv"

NPIX = 640000
# each pix: 2 pc x 2 pc
# Initial cloud size pi * 100 pc x 100 pc ~ 7854 pix
NPIXCLOUD = 7854

if(PROJ == "perpendicular"):
    NPIXLIM = (int)(NPIXCLOUD * 1.40) # x300v1700c
if(PROJ == "parallel"):
    NPIXLIM = (int)(NPIXCLOUD * 0.16) # x300v1700c

# Use 2 * NPIXCLOUD;
# Sort column densities from large to small, until 2 x NPIXCLOUD are found

# 1   2   3  4   5      6    7    8     9    10
# OVI CIV NV CII NeVIII CIII MgII SiIII SiIV HI
ionnames = ["XX","OVI","CIV","NV","CII","NeVIII","CIII","MgII","SiIII","SiIV","HI"]
ionorder = [-1, 5, 3, -1, -1, 6, 2, 7, -1, 8, 0]

ionidxs = [10, -1, 6, 2, -1, 1, 5, 7, 9] # Mapping from the 9ions to ionnames
clrs = ["blue", "steelblue", "lime", "green", "plum", "orange", "red", "purple", "olive"]
lgds = ["HI","HeII(N/A)","CIII","CIV","OIV(N/A)","OVI","NeVIII","MgII","SiIV"]
ions = ["HI","HeII","CIII","CIV","OIV","OVI","NeVIII","MgII","SiIV"]

# 0  1    2    3   4   5   6      7    8
# HI HeII CIII CIV OIV OVI NeVIII MgII SiIV
# No: HeII, OIV; Extra: SiIII, NV, CII

# Format:
# Row: 9 ions
# column: PDF for each ion, 100 cells (101 lines)
NCELLS = 100
cols = range(1, 11)

def generate_pdf_tables():
    tabout = []
    if(PROJ == "perpendicular"):
        ioncds = ioformat.rcol(fnamex, cols, separator=",", linestart=1)
        foutname = "pdfs/"+modelstr+"_i9_x_pdf.dat"
    else:
        ioncds = ioformat.rcol(fnamey, cols, separator=",", linestart=1)
        foutname = "pdfs/"+modelstr+"_i9_y_pdf.dat"
    edges = linspace(0., 1., NCELLS+1)
    # pdb.set_trace()
    for i in range(9):
        idx = ionidxs[i] 
        if(idx == -1):
            tabout.append(array([0.0] * (NCELLS+1)))
            print ("---")
            continue
        ionname = ionnames[idx]
        print (ionname)
        ioncds[idx-1].sort()
        prob = log10(ioncds[idx-1][::-1][:NPIXLIM][::NPIXLIM/NCELLS][:NCELLS+1])
        tabout.append(prob)

    fout = open(foutname, "w")
    fout.write("#P(>logN) HI HeII CIII CIV OIV OVI NeVIII MgII SiIV\n")
    for i in range(NCELLS+1):
        line = "%5.3f " % (edges[i])
        for j in range(9):
            line += "%6.3f " % (tabout[j][i])
        line += "\n"
        fout.write(line)
    fout.close()

def draw():
    fig, ax = plt.subplots(1,1,figsize=(8,6))
    # models = ['T0.3_v1700_chi300_cond_0-01',\
    #           'T0.3_v1700_chi300_cond_0-1',\
    #           'x300v1700c',\              
    #           'T0.3_v1700_chi300_cond_10',\
    #           'T0.3_v1700_chi300_cond_100']
    # lgds_uv = ['x0.01', 'x0.1', 'x1', 'x10', 'x100']
    models = ['T0.3_v1700_chi300_cond_10',\
              'T0.3_v1700_chi300_cond_100',\
              'T0.3_v1700_chi300_cond_1e5',\
              'T0.3_v1700_chi300_cond_1e6',\              
              'T0.3_v1700_chi300_cond_1e7']
    lgds_uv = ['x10', 'x100', 'x1e5', 'x1e6', 'x1e7']
    # lstyles = [":", "--", "-", "--", ":"]
    lstyles = ["dotted", "dashed", (0,(5,1)), (0,(3,1,1,1)), "-"]    
    lw = [1,1,2,1,1]
    foutnames = []
    for model in models: foutnames.append("pdfs/"+model+"_i9_y_pdf.dat")
    # Ncorr = 0.16 / 1.4
    Ncorr = 1.0
    ilist = range(9)
    # ilist = [0,3,5,7]
    for fi, foutname in enumerate(foutnames):
        print (foutname)
        tab = genfromtxt(foutname, names=True)
        edges = tab['PlogN']
        for i in ilist:
            ion = ions[i]
            ax.plot(tab[ion], edges, linestyle=lstyles[fi], color=clrs[i], linewidth=lw[fi])
            # ax.plot(taby[ion] + log10(Ncorr), edges, "--", color=clrs[i])
    ax.set_xlim(12., 22.)
    ax.set_ylim(0., 1.)
    ax.set_xlabel(r"$\log(N_{ion})$")
    ax.set_ylabel(r"P")
    ax.set_title("x300v1700")
    # if(PROJ == "perpendicular"):
    #     ax.set_title("x300v1700-lc, Perpendicular")
    # else:
    #     ax.set_title("x300v1700-lc, Parallel")
    from pltastro import legend
    lgd = legend.legend(ax)
    for i in ilist:
        lgd.addLine((lgds[i], clrs[i], "-", 1))
    lgd.loc = "lower left"
    lgd.draw()
    lgd2 = legend.legend(ax)
    for i in range(5):
        lgd2.addLine((lgds_uv[i], "black", lstyles[i], lw[i]))
    lgd2.loc = "upper right"
    lgd2.draw()
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()

print ("done")
