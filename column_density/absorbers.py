from numpy import genfromtxt, log10
import matplotlib.pyplot as plt
import matplotlib as mpl

fbase = "/scratch/shuiyao/los/l25n288-phew-m5/"
angles = range(10, 20, 2)
fig, ax = plt.subplots(1,1,figsize=(7,6))
cmap = plt.get_cmap('jet')
for angle in angles:
    f = fbase+"cloudszfile.l25n288-phew-m5."+str(angle)+".50_0"
    tab = genfromtxt(f, names=("IonID","logN","Z","NCorr","Mc","Rc"))
    tab = tab[tab['IonID'] == 0.0]
    clrs = cmap(tab['Mc'])
    ax.scatter(tab['logN']+log10(tab['NCorr']), tab['Rc'], c=clrs, s=6)

ax.yaxis.set_ticks([0.1, 1.0, 10.0])
ax.yaxis.set_ticklabels("0.1", "1.0", "10.0")
ax.set_ylim(0.1, 10.0)
plt.xlabel("logN(HI)")
plt.ylabel("Rc [kpc]")
plt.yscale("log")

axcbar = fig.add_axes([0.15,0.1,0.7,0.015])
norm1 = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.cm.jet, norm=norm1, orientation="horizontal")
cdcbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
cdcbar.set_ticklabels(["0.0","0.2","0.4","0.6","0.8", "1.0"])
cdcbar.set_label("Mc/Mc_0")

plt.subplots_adjust(bottom=0.25)

plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.close()

print "Done."
