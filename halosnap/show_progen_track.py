import ioformat
import matplotlib.pyplot as plt

fname = "progen_02253"
x, y, z = ioformat.rcol(fname, [1,2,3])
plt.plot(x, z, "b.-")
for i in range(len(x)):
    plt.text(x[i], y[i], str(i), fontsize=6)
plt.axis([-0.5, 0.5, -0.5, 0.5])
plt.show()

