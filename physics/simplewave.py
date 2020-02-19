# Formulae for a simple wave with a polytropic index n_poly
# We would parameterize everything in terms of x = 1 - u/u_esc
# u_esc is a function of n_poly only, so is x

# Temperature profile depends on n_poly, but density/pressure are very robust

# rho(x) = rho_0 * x ** (2/(n-1))
# P(x) = P_0 * x ** (2n/(n-1))
# x = 1 - u/u_esc is actually c / c_0

from numpy import inf, exp, array

def x_par(u, cs, n_poly=5./3.):
    if(n_poly > 1.0):
        b = 2.0 / (n_poly - 1.0)
        return 1 - (u / cs) / b
    elif(n_poly == 1.0): # a different diffinition
        return u / cs
    else:
        print "Error: n < 1.0."        

def fac_density(x, n_poly=5./3.):
    if(n_poly > 1.0):    
        b = 2.0 / (n_poly - 1.0)
        return x ** b
    elif(n_poly == 1.0):
        return exp(-x)
    else: print "Error: n < 1.0."        

def fac_pressure(x, n_poly=5./3.):
    if(n_poly > 1.0):    
        b = 2.0 / (n_poly - 1.0)
        return x ** (n_poly * b)
    elif(n_poly == 1.0):
        return exp(-x)
    else: print "Error: n < 1.0."        

def u_esc(cs, n_poly=5./3., isocut=5.0):
    if(n_poly > 1.0):    
        return 2.0 * cs / (n_poly - 1.0)
    elif(n_poly == 1.0):
        return isocut * cs # an arbitrary value
    else: print "Error: n < 1.0."        

def show(n_poly=5./3., color1="blue", color2="red"):
    from numpy import linspace
    import matplotlib.pyplot as plt
    u = linspace(0., u_esc(1., n_poly), 100)
    x = x_par(u, 1., n_poly)
    rho = fac_density(x, n_poly)
    P = fac_pressure(x, n_poly)
    if(n_poly > 1): T = x ** 2
    else: T = array([1.0] * len(x))    
    plt.plot(u, rho, "-", color=color1)
    plt.plot(u, T, "-", color=color2)
    plt.plot(u, P, "--", color=color1)

def check():
    import matplotlib.pyplot as plt
    cmap1 = plt.get_cmap("Blues")
    cmap2 = plt.get_cmap("Reds")
    ns = [1.0, 1.2, 1.4, 5./3.]
    clri = [50, 100, 150, 200]
    for i in range(4):
        show(ns[i], color1=cmap1(clri[i]), color2=cmap2(clri[i]))
    plt.axis([0.0, 10.0, 0.0, 1.0])
    plt.show()
    
