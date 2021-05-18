import ioformat
a = ioformat.rcol("./outputs_z108.txt", [0])

fout = open("redshifts.txt", "w")
fout.write("#snapnum a zred\n")
for i in range(len(a)):
    z = 1./a[i] - 1.
    line = "%3d %7.5f %7.5f\n" % (i, a[i], z)
    fout.write(line)
fout.close()

print ("DONE.")
