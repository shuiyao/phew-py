from myinit import *

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 30
redges = linspace(0., 1., nbins+1)
dr = 1. / nbins

def find_radial_bin(r):
    if(r > 1.0): return -1
    bin_idx = r / dr
    return (int)(bin_idx)

REDSHIFT = 0.25
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 0.25): zstr = "098"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"

def calculate_gas_profile(model, zstr):
    for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
        # fname << 
        fname = DIRS['SCIDATA'] + model + "/snapshot_"+zstr+".gas."+mstr
        foutname = DIRS['SCIDATA'] + model + "/gasprof_"+zstr+"_"+mstr

        print ("Reading: ", fname)
        print ("Generating New Mass Profile File.")
        tab = genfromtxt(fname, names=True)

        mass = array([0.0] * (nbins + 1))
        mwind = array([0.0] * (nbins + 1))
        mmix = array([0.0] * (nbins + 1))
        mwcold = array([0.0] * (nbins + 1))
        mwhot = array([0.0] * (nbins + 1))                
        mcold = array([0.0] * (nbins + 1))
        mhot = array([0.0] * (nbins + 1))
        mism = array([0.0] * (nbins + 1))
        mmaxhot = array([0.0] * (nbins + 1))
        mmaxcold = array([0.0] * (nbins + 1))
        mmaxwhot = array([0.0] * (nbins + 1))
        mmaxwcold = array([0.0] * (nbins + 1))
        mzcold = array([0.0] * (nbins + 1))
        mzhot = array([0.0] * (nbins + 1))
        mzism = array([0.0] * (nbins + 1))
        mzwind = array([0.0] * (nbins + 1))        

        for i in range(len(tab)):
            part = tab[i]
            bidx = find_radial_bin(part['dr'] / part['Rvir'])
            mass[bidx] += part['Mass'] # Total Mass
            mmix[bidx] += part['WMass'] # Total Mixed Mass
            if(part['SfFlag'] == 1):
                mism[bidx] += part['Mass']
                mzism[bidx] += part['Mass'] * part['Z']
            else:
                if(part['Mc'] > 0): # A PhEW Particle
                    mwind[bidx] += part['Mass']
                    mmix[bidx] -= part['WMass']
                    mzwind[bidx] += part['Mass'] * part['Z']                    
                if(part['Mc'] == 0): # A normal gas particle
                    if(part['logT'] > 5.5):
                        mhot[bidx] += part['Mass'] - part['WMass']
                        mwhot[bidx] += part['WMass']
                        mzhot[bidx] += part['Mass'] * part['Z']
                    else:
                        mcold[bidx] += part['Mass'] - part['WMass']
                        mwcold[bidx] += part['WMass']
                        mzcold[bidx] += part['Mass'] * part['Z']                        
                    if(part['Tmax'] > 5.5):
                        mmaxhot[bidx] += part['Mass'] - part['WMass']
                        mmaxwhot[bidx] += part['WMass']
                    else:
                        mmaxcold[bidx] += part['Mass'] - part['WMass']
                        mmaxwcold[bidx] += part['WMass']
                if(part['Mc'] < 0): # Tricky, now very rare
                    mwind[bidx] += part['Mass'] - part['WMass']
        fwind = mwind/mass
        fout = open(foutname, "w")
        fout.write("#r Mass Mcold Mhot Mwind Mwcold Mwhot Mism Mmaxcold Mmaxhot Mmaxwcold Mmaxwhot Mzcold Mzhot Mzism Mzwind\n")
        for i in range(len(redges)):
            line = "%5.3f %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n" % \
                (redges[i], \
                 mass[i], mcold[i], mhot[i], \
                 mwind[i], mwcold[i], mwhot[i], mism[i], \
                 mmaxcold[i], mmaxhot[i], mmaxwcold[i], mmaxwhot[i], \
                 mzcold[i], mzhot[i], mzism[i], mzwind[i])
            fout.write(line)
        fout.close()

calculate_gas_profile("l50n288-phew-m5-spl", zstr) 
print ("Done.")
