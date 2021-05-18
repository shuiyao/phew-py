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
import os

import pandas as pd

print ("compiled.")
SYSTEM = 'unity'
ATOL = 1.e-4 # Tolerance for ascale (a(SF) < a(Acc) + ATOL)
# DTTOL = 1000. * ac.myr
DTTOL = 1.e9 # yr

STELLAR_AGE_IN_ASCALE = True # Format of the P(Star).Age
SFRINFO_FORMATS = ["Gadget3", "GIZMO-PhEWOff", "GIZMO-PhEW-Extra"]

import sys
modelname = sys.argv[1]
snapstr = sys.argv[2]
lbox = (float)(sys.argv[3])
NCPU = (int)(sys.argv[4])
flag = (int)(sys.argv[5])
if(len(sys.argv) > 6): mstr = "."+sys.argv[6]
else: mstr = ""

unit_m = 3469578.81574 * (lbox / 50.) ** 3
print ("unit_m =", unit_m)

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

def format_of_sfrinfo_file(flag_format):
    if(flag_format == "GIZMO-PhEW-Extra"): # update: 20200724
        fformat={'names': ('atime', 'AccKey', 'alast', 'tmax', 'mass', 'wmass', 'Z', 't1', 't3'),
                 'formats': ('f8', 'i8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}
        cols = (0,2, 3,4,5,6,7, 8,10)
    if(flag_format == "GIZMO-PhEWOff"): # update: 20210107
        fformat={'names': ('atime', 'AccKey', 'alast', 'tmax', 'mass', 'Z'),
                 'formats': ('f8', 'i8', 'f8', 'f8', 'f8', 'f8')}
        cols = (0,2,3,4,5,7)
    if(flag_format == "Gadget3"):
        fformat={'names': ('atime', 'ID', 'alast', 'tmax', 'mstar', 'mass', 'Z'),
                 'formats': ('f8', 'i8', 'f8', 'f8', 'f8', 'f8', 'f8')}
        cols = (0,1,2,3,5,6,7)
    return fformat, cols

def read_stars_data(fname, stellar_age_in_ascale = True):
    stars = genfromtxt(fname, dtype='i8,i8,i8,i8,f8,f8,f8', names=True)
    # Idx, ID, GID, HID, Mass, Tmax, Age
    print (stars.dtype)
    stars['Mass'] = stars['Mass'] * unit_m * 1.e10 / 0.7 * 2.0e33
    # 'Age' used to be the real age. Now is the a_form
    if(stellar_age_in_ascale == False): # Force Convert the 'Age' field to ascale
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

def load_sfrinfo(stars, gals, sfrinfobase, outname, flag_format):
    # flag_format: The format of the SFRINFO file
    # Let's be concerned ONLY with STARS.
    # Idea: Loop over all sfrinfo.* files.
    # For each sfrinfo.* file, guided loop over particle ID.
    #   - a(sfrinfo) < a(fstars) is a sure thing.
    # Both files should be ranked with PID.
    #   - For a fast searching, should bin PID into bins
    # Each time sfrinfo.* finds a likely match for a PID, find
    #   - tmax(sfrinfo), tmax(fstars)
    # If tmax match, then YES!

    fformat, cols = format_of_sfrinfo_file(flag_format)
    print ("flag_format = ", flag_format)
    print (fformat)
    print (cols)

    # stars = sort(stars, order='ID') # Order By StarP.ID(AccKey)
    acckey_to_sidx = dict()
    acckey_to_alast = dict()
    acckeys = stars['ID']
    smass = stars['Mass']
    for i in range(len(acckeys)):
        acckey_to_sidx[acckeys[i]] = i
        acckey_to_alast[acckeys[i]] = 0.0
    
    nparts = len(stars)
    aacc, alast, pmass, wmass, metals = [-1]*nparts, [-1]*nparts, [0.0]*nparts, [0.0]*nparts, [-1.0]*nparts # update: 20200110
    if(flag_format == "GIZMO-PhEW-Extra"): tavg, sigavg = [-1]*nparts, [0.0]*nparts
    # ---------------- LOOP: NCPU ----------------
    for fi in range(NCPU):
        fsfrinfo = sfrinfobase + "sfrinfo." + str(fi)
        print ("Reading: ", fsfrinfo)
        if(not os.path.exists(fsfrinfo)): break
        # acc = genfromtxt(fsfrinfo, dtype=fformat, usecols=cols) # update: 20200110
        acc = pd.read_csv(fsfrinfo, usecols=['atime', 'AccKey', 'LastSFTime', 'Tmax', 'WindMass', 'Mass', 'Z', 't1', 't3'])
        acc.rename(columns={"LastSFTime":"alast", "Tmax":"tmax", "WindMass":"wmass", "Mass":"mass"}, inplace=True)
        acc = acc[acc['alast'] == 0] # excluding 'spurious accretions'
        acckeys = acc['AccKey']
        for iacc in range(len(acc)):
            last_sf_time = acc['alast'].iloc[iacc] # LastSFTime
            thiskey = acckeys.iloc[iacc]
            if(thiskey in acckey_to_sidx): # Match
                # Now find the most recent major accretion event
                if(abs(acc['alast'].iloc[iacc]) >= acckey_to_alast[thiskey]):
                    acckey_to_alast[thiskey] = abs(acc['alast'].iloc[iacc])
                else: continue
                sidx = acckey_to_sidx[acckeys.iloc[iacc]]
                aacc[sidx] = acc['atime'].iloc[iacc] # Accretion Time
                alast[sidx] = last_sf_time
                pmass[sidx] = acc['mass'].iloc[iacc] # Particle Mass At Accretion
                metals[sidx] = acc['Z'].iloc[iacc] # Metallicity
                if(flag_format == "GIZMO-PhEWOff"):
                    wmass[sidx] = acc['mass'].iloc[iacc] if last_sf_time < 0 else 0.0
                if(flag_format == "GIZMO-PhEW-Extra"):
                    wmass[sidx] = acc['wmass'].iloc[iacc] # add: 20200110
                    tavg[sidx] = acc['t1'].iloc[iacc] / acc['wmass'].iloc[iacc] # add: 20200724
                    if(tavg[sidx] > 1.34e10):
                        tavg[sidx] = 1.0
                    elif(tavg[sidx] < 1.e6):
                        tavg[sidx] = 0.0
                    else:
                        tavg[sidx] = acosmic(tavg[sidx])
                    sigavg[sidx] = acc['t3'].iloc[iacc] / acc['wmass'].iloc[iacc] # add: 20200724
    fout = open(outname, "w")
    if(flag_format == "GIZMO-PhEW-Extra"):
        fout.write("#a_form a_acc a_last Mass StarMass WindMass WindAge WindSig Tmax Z GID HID\n") # update: 20200724
    else:
        fout.write("#a_form a_acc a_last Mass StarMass WindMass Tmax Z GID HID\n") # update: 20200110        
    for istars in range(len(stars)):
        gid, hid = stars['GID'][istars], stars['HID'][istars]
        # update: 20200110
        if(flag_format == "GIZMO-PhEW-Extra"):        
            line = "%7.5f %7.5f % 7.5f %7.5e %7.5e %7.5e %7.5f %5.1f % 5.3f %7.5e %5d %5d\n" % \
            (stars['Age'][istars], aacc[istars], alast[istars], \
             pmass[istars], smass[istars], wmass[istars], tavg[istars], sigavg[istars], \
             stars['Tmax'][istars], metals[istars], \
             gid, hid)
        else:
            line = "%7.5f %7.5f % 7.5f %7.5e %7.5e %7.5e % 5.3f %7.5e %5d %5d\n" % \
            (stars['Age'][istars], aacc[istars], alast[istars], \
             pmass[istars], smass[istars], wmass[istars], \
             stars['Tmax'][istars], metals[istars], \
             gid, hid)
        fout.write(line)
    fout.close()

if(SYSTEM == "eagle"):
    galbase = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"
    fname = galbase+modelname+"_"+snapstr+".stars"+mstr
    outname = galbase+modelname+"_"+snapstr+".starinfo"+mstr
    sfrinfobase = "/proj/shuiyao/"+modelname+"/SFRINFO/"
    skidbase = "/proj/shuiyao/"+modelname+"/"
if(SYSTEM == "unity"):
    galbase = "/home/shuiyao_umass_edu/scidata/"+modelname+"/"
    fname = galbase+modelname+"_"+snapstr+".stars"+mstr
    outname = galbase+modelname+"_"+snapstr+".starinfo"+mstr
    sfrinfobase = "/nas/astro-th/shuiyao/"+modelname+"/SFRINFO/"
    skidbase = "/nas/astro-th/shuiyao/"+modelname+"/"
stars = read_stars_data(fname, STELLAR_AGE_IN_ASCALE)
gals = read_gal_data(skidbase, snapstr)
load_sfrinfo(stars, gals, sfrinfobase, outname, SFRINFO_FORMATS[flag])

