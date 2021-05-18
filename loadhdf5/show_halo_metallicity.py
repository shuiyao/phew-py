import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from myinit import *

#model = "l50n288-phew-m5"
models = ["l50n288-phewoff", "l50n288-phew-m5", "l50n576-phew-m5"]
lstyles = [':', '--', '-']

zstr = "098"

fig, axs = plt.subplots(3,1,figsize=(6,10))

def plot_model(model, axs, lstyle="-"):
    f = os.path.join(DIRS['SCIDATA'], model, "Zhist_"+zstr)
    print("plotting: ", f)
    df = pd.read_csv(f)
    df = df.drop(index=0)
    axs[0].plot(df.logZ, df.mh11c, color="blue", linestyle=lstyle)
    axs[0].plot(df.logZ, df.mh11h, color="red", linestyle=lstyle)
    axs[0].plot(df.logZ, df.mh11c+df.mh11h, color="black", linestyle=lstyle)    
    axs[1].plot(df.logZ, df.mh12c, color="blue", linestyle=lstyle)
    axs[1].plot(df.logZ, df.mh12h, color="red", linestyle=lstyle)
    axs[1].plot(df.logZ, df.mh12c+df.mh12h, color="black", linestyle=lstyle)        
    axs[2].plot(df.logZ, df.mh13c, color="blue", linestyle=lstyle)
    axs[2].plot(df.logZ, df.mh13h, color="red", linestyle=lstyle)
    axs[2].plot(df.logZ, df.mh13c+df.mh13h, color="black", linestyle=lstyle)        

for i, model in enumerate(models):
    plot_model(model, axs, lstyle=lstyles[i])

for i in range(3):
    axs[i].set_xlim(-4.0, 1.0)
    axs[i].set_ylabel("f(Z)")
    axs[i].set_ylim(0.0, 0.10)
    axs[i].set_yticks([0.0, 0.03, 0.06, 0.09])
axs[0].set_xticks([])
axs[1].set_xticks([])
axs[2].set_xlabel(r"$\log(Z/Z_\odot)$")
axs[2].set_xticks([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0])
axs[0].text(0.05, 0.85, r"$11.0 < \log(M_{vir}/M_\odot) < 11.5$", transform=axs[0].transAxes)
axs[1].text(0.05, 0.85, r"$11.85 < \log(M_{vir}/M_\odot) < 12.15$", transform=axs[1].transAxes)
axs[2].text(0.05, 0.85, r"$12.85 < \log(M_{vir}/M_\odot) < 13.15$", transform=axs[2].transAxes)
axs[0].set_title('z = 0.25')
fig.subplots_adjust(hspace=0.01, left=0.18)

from pltastro import legend
lgds = []
lgd = legend.legend(axs[0])
lgd.addLine(("cold", "blue", '-', 1))
lgd.addLine(("hot", "red", '-', 1))
lgd.addLine(("all", "black", '-', 1))
lgd.loc = 'upper right'
lgd.draw()
lgd2 = legend.legend(axs[2])
for i, model in enumerate(models):
    lgd2.addLine((model, "black", lstyles[i], 1))
lgd2.loc = 'center left'
lgd2.fontsize = 10
lgd2.draw()

plt.savefig(os.path.join(DIRS['FIGURE'], 'tmp.pdf'))

plt.show()
