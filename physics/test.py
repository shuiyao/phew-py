import matplotlib.pyplot as plt

print "compiled."
fig = plt.figure(0, figsize=(7,3))
nmax = 40
c = range(1,nmax)
y = []
count = 0
for ci in c:
    b = range(ci + 1, nmax)
    for bi in b:
        a = range(bi + 1, nmax)
        for ai in a:
            y.append(ai * bi * ci * (100 * ai + 10 * bi + ci))
            if(y[-1] < 500000): count += 1
            plt.plot([y[-1], y[-1]], [0,1], "b-")

# plt.plot(y, y, "b.")
# plt.yscale("log")
print count
plt.axis([1,500000,0,1])
# plt.xscale("log")
plt.show()
# plt.hist(y, bins=100)            
print "finished."            
        
