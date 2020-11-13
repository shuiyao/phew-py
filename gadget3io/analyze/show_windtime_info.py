from mymod import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic

# Show wind time information by redshifts z = 0.25, z = 1.0, z = 2.0

hparam = 0.7
unit_m = 1.e10 / hparam

CUMULATIVE_DISTRIBUTION = True

models = ["l25n144-phew-m5-spl", "l25n288-phew-m5-spl"]
lgds = ["25/144,Split", "25/288,Split"]

lstyles = ["--", "-"]

from pltastro import frame, draw
import config_mpl
frm = frame.multi(3,1)
pars = frm.params
pars.figsize = (5, 9)
pars.left = 0.2
pars.top = 0.92
pars.bottom = 0.2
panels = frm.panels
panels.set_xlabels(r"$Age [Gyr]$")
panels.set_ylabels("Fraction")
# panels.set_xlims(0.0, 1.0)
panels.set_xlims(0.0, 10.0)
if(CUMULATIVE_DISTRIBUTION):
    panels.set_ylims(0.0, 1.0)
else:
    panels.set_ylims(0.0, 0.4)
    panels.set_yticks([0.0, 0.1, 0.2, 0.3])

fig, axs = draw(frm)

captions = [
    r"$11.0 < M_\mathrm{vir} < 11.5$",\
    r"$11.85 < M_\mathrm{vir} < 12.15$",\
    r"$12.85 < M_\mathrm{vir} < 13.15$"\    
]
mmin = [11.0, 11.85, 12.85]
mmax = [11.5, 12.15, 13.15]

clrs = ["blue", "orange", "red", "grey"]
cols = ["Mh11", "Mh12", "Mh13", "All"]
zred = REDSHIFT

def plot_for_redshift(axi, zstr, zred):
    for modeli in range(len(models)):
        fname = "/home/shuiyao_umass_edu/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".windage"
        tab = genfromtxt(fname, names=True)
        age = (tcosmic(1./(zred+1.)) - tcosmic(tab['atime'])) / 1.e9
        for i in range(4):
            frac = cumsum(tab[cols[i]])
            frac = frac / frac[-1]
            axs[axi].plot(age, frac, color=clrs[i], linestyle=lstyles[modeli])
        
axs[0].set_title("Wind Time, z="+str(REDSHIFT)[:4])

# plot_for_redshift(0, "108", 0.0)
plot_for_redshift(1, "078", 1.0)

from pltastro import legend
lgd1 = legend.legend(axs[0])
lgd1.loc = "lower right"
lgd1.addLine((lgds[0], "black", lstyles[0], 1))
lgd1.addLine((lgds[1], "black", lstyles[1], 1))
lgd1.draw()
# lgd2 = legend.legend(axs[2])
# lgd2.loc = "lower right"
# lgd2.addLine(("cold wind", "cyan", "-", 1))
# lgd2.addLine(("hot wind", "magenta", "-", 1))
# lgd2.addLine(("total wind", "orange", "-", 1))
# lgd2.draw()

plt.savefig("/home/shuiyao_umass_edu/figures/tmp.pdf")
plt.show()

print ("Done.")
