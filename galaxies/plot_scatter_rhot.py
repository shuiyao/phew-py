import matplotlib.pyplot as plt
import ioformat

fbase = "/proj/shuiyao/"
model = "l25n144-phew"
fname = fbase + model + "/snapshot_108"

# rho = ioformat.rcol(fname+".rho", [0])
# temp = ioformat.rcol(fname+".temp", [0])
# delayt = ioformat.rcol(fname+".delayt", [0])

# rhow, tempw = [], []
# for i in range(rho.__len__()):
#     if(delayt[i] > 0):
#         rhow.append(rho[i])
#         tempw.append(temp[i])

plt.plot(rho[::50], temp[::50], "b.", markersize=3, alpha=0.2)
plt.plot(rhow, tempw, "r.", alpha=0.2, markersize=1)
plt.xscale("log")
plt.yscale("log")
plt.show()

