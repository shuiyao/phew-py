# Compile Information for star particles.
# Including their:
# Mstar, Mvir, aform, aacc, Tmax, alast

# Sfrinfo:
# a, ID, LastSFTime, Tmax, ...
# LastSFTime =
#   - 0: First Time SF, pristine gas accretion
#   - +: non-SF -> SF
#   - -: Wind -> SF
# Tmax is set to 0 at launch.
from cosmology import acosmic, tcosmic
import ioformat
from numpy import genfromtxt, linspace, array, sort, insert
from numpy import where, inf, log10
import config_mpl
import matplotlib.pyplot as plt
import config_mpl
from astroconst import pc, ac

print "compiled."
unit_m = 3469578.81574
NCPU = 256
ATOL = 1.e-4 # Tolerance for ascale (a(SF) < a(Acc) + ATOL)
# The criteria for discarding spurious accretion:
# If dM* < DMTOL: It's very likely spurious
#   - However, if in this case dt > 1 Gyr, we still take them into account
# If dM* > DMTOL: We will base solely on dt < 200 Myr to discard them
DTTOL = 1000. * ac.myr
DMTOL = 0.2
NEGATIVE_TMAX_FOR_WINDS = False
# Note that in some later version of the Gadget code (e.g., p50n288dsw), I make the Tmax of particles negative if their last accretion is from winds. But in previous versions it's not there.
modelname = "p50n288sw"
snapstr = "058"

class skidgal():
    def __init__(self, mstar, mvir):
        self.mvir = log10(mvir / 0.7)
        self.mstar = log10(mstar * unit_m * 1.e10 / 0.7)
        
# read_gal_data(fname)
# show_sfh(gal, fraction=False):

# load_sfrinfo(): The primary function to combine the .stars and the SFRINFO data

def read_gal_data(skidbase, snapstr):
    fgal = "%s/gal_z%s.stat" % (skidbase, snapstr)
    fso = "%s/so_z%s.sovcirc" % (skidbase, snapstr)
    mstar = ioformat.rcol(fgal, [4])
    mvir = ioformat.rcol(fso, [1], linestart=1)
    gals = []
    for i in range(len(mstar)):
        gals.append(skidgal(mstar[i], mvir[i]))
    return gals

def read_stars_data(fname):
    stars = genfromtxt(fname, dtype='i8,i8,i8,i8,f8,f8,f8', names=True)
    # Idx, ID, GID, HID, Mass, Tmax, Age
    print stars.dtype
    stars['Mass'] = stars['Mass'] * unit_m * 1.e10 / 0.7
    for i in range(len(stars['Age'])):
        if(stars['Age'][i] > 0):
            stars['Age'][i] = acosmic(stars['Age'][i])
    return stars

def read_galacc_data(fname): # The compiled file.
    galacc = genfromtxt(fname, dtype='i8,f8,f8,f8,f8,f8,f8', names=True)
    # Mass in unit of Msolar
    # "#PID a_form a_acc a_last Mass Tmax\n"    
    return galacc

def match_check(a, b, tol=0.2): # Strange the tmax does not always match exactly
    if(abs(abs(a) - abs(b)) < tol): return True
    else: return False

def discard_spurious_accretion(atime, alast):
    # ac is THIS accretion event to be evaluated
    if(alast <= 0): return 0
    if(tcosmic(atime) - tcosmic(alast) < DTTOL): return 1
    return 0

def load_sfrinfo(stars, gals, sfrinfobase, outname):
    # Let's be concerned ONLY with STARS.
    # Idea: Loop over all sfrinfo.* files.
    # For each sfrinfo.* file, guided loop over particle ID.
    #   - a(sfrinfo) < a(fstars) is a sure thing.
    # Both files should be ranked with PID.
    #   - For a fast searching, should bin PID into bins
    # Each time sfrinfo.* finds a likely match for a PID, find
    #   - tmax(sfrinfo), tmax(fstars)
    # If tmax match, then YES!
    fformat={'names': ('atime', 'ID', 'alast', 'tmax', 'mstar'),
             'formats': ('f8', 'i8', 'f8', 'f8', 'f8')}
    # fsfrinfo = sfrinfobase + "sfrinfo.56"
    stars = sort(stars, order='ID')
    nparts = len(stars)
    aacc, alast, tmaxdiff, ms = [-1]*nparts, [-1]*nparts, [10.]*nparts, [-1]*nparts
    for fi in range(NCPU):
        fsfrinfo = sfrinfobase + "sfrinfo." + str(fi)
        print "Reading: ", fsfrinfo
        acc = genfromtxt(fsfrinfo, dtype=fformat, usecols=(0,1,2,3,5))
        acc = sort(acc, order='ID')
        istars, nextPID = 0, stars['ID'][0]
        for iacc in range(len(acc)): # Loop: SFRINFO
            if(istars > nparts - 1): break
            if(acc['ID'][iacc] < nextPID): continue
            # Now acc.ID >= nextPID
            # while(stars['ID'][istars] == nextPID): # Start counting in the stars array
            while(acc['ID'][iacc] >= nextPID): # Start counting in the stars array  
                if(acc['ID'][iacc] == nextPID): # For each stars particle, see if a match is found
                    # Do the spurious accretion check
                    if(not discard_spurious_accretion(acc['atime'][iacc], acc['alast'][iacc])):
                        if(match_check(acc['tmax'][iacc], stars['Tmax'][istars], tol=tmaxdiff[istars])): # tmaxs match. When Tmax matches, it restrict the star to its closest major accretion event, whether it's primordial or wind
                            if(acc['atime'][iacc] < stars['Age'][istars]+ATOL): # SF after accretion
                                aacc[istars] = acc['atime'][iacc]
                                alast[istars] = acc['alast'][iacc]
                                ms[istars] = acc['mstar'][iacc]         
                                tmaxdiff[istars] = abs(abs(acc['tmax'][iacc])-abs(stars['Tmax'][istars]))
                        # End (acc['ID'][iacc] == nextPID)
                istars += 1
                if(istars > nparts - 1): break # No particle to match!            
                # while(acc['ID'][iacc] >= nextPID):
                nextPID = stars['ID'][istars] # now updated.
    # Now we got aacc, alast, tmaxdiff.
    # Now Write:
    fout = open(outname, "w")
    fout.write("#a_form a_acc a_last Mass Mstar Tmax GID HID\n")
    for istars in range(len(stars)):
        if(NEGATIVE_TMAX_FOR_WINDS):
            if(alast[istars] < 0): stars['Tmax'][istars] *= -1
        gid, hid = stars['GID'][istars], stars['HID'][istars]
        if(gid != 0): mstar = gals[gid-1].mstar
        else: mstar = -inf
        if(hid != 0): mvir = gals[abs(hid)-1].mvir
        else: mvir = -inf
        if(gid != hid): # non-central
            mvir *= -1
        line = "%7.5f %7.5f % 7.5f %7.5e %5.2f % 5.3f %5d %5d\n" % \
        (stars['Age'][istars], aacc[istars], alast[istars], \
         stars['Mass'][istars], ms[istars], stars['Tmax'][istars],
         gid, hid)
        fout.write(line)
    fout.close()

import sys
modelname = sys.argv[1]
snapstr = sys.argv[2]
mstr = "."+sys.argv[3]
if(modelname == "p50n288dsw"): NEGATIVE_TMAX_FOR_WINDS = True
galbase = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"
fname = galbase+modelname+"_"+snapstr+".stars"+mstr
outname = galbase+modelname+"_"+snapstr+".starinfo"+mstr
sfrinfobase = "/scratch/shuiyao/data/"+modelname+"/SFRINFO/"
skidbase = "/scratch/shuiyao/data/"+modelname+"/"
stars = read_stars_data(fname)
gals = read_gal_data(skidbase, snapstr)
load_sfrinfo(stars, gals, sfrinfobase, outname)
