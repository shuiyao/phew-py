# Analyze the history of one galaxy.

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
from numpy import where
import config_mpl
import matplotlib.pyplot as plt
import config_mpl

print "compiled."
unit_m = 3469578.81574
NCPU = 256
ATOL = 1.e-4 # Tolerance for ascale (a(SF) < a(Acc) + ATOL)
fname = "p50n288sw_gals.lst"
fbase = "/scratch/shuiyao/data/p50n288sw/SFRINFO/"
galbase = "/scratch/shuiyao/sci/gadget3io/p50n288sw/"
fgallst = "p50n288zw_p50n288sw.lst"

# read_gal_data(fname)
# show_sfh(gal, fraction=False):
# load_sfrinfo():

def read_gal_data(fname):
    gal = genfromtxt(fname, dtype='i8,i8,f8,f8,f8', names=True)
    # Idx, ID, Mass, Tmax, Age
    print gal.dtype
    gal['Mass'] = gal['Mass'] * unit_m * 1.e10 / 0.7
    for i in range(len(gal['Age'])):
        if(gal['Age'][i] > 0):
            gal['Age'][i] = acosmic(gal['Age'][i])
    return gal

def read_galacc_data(fname): # The compiled file.
    galacc = genfromtxt(fname, dtype='i8,f8,f8,f8,f8,f8,f8', names=True)
    # Mass in unit of Msolar
    # "#PID a_form a_acc a_last Mass Tmax\n"    
    return galacc

def match_check(a, b, tol=0.2): # Strange the tmax does not always match exactly
    if(abs(a - b) < tol): return True
    else: return False

def load_sfrinfo(gal, sfrinfobase, outname):
    # Let's be concerned ONLY with STARS.
    # Idea: Loop over all sfrinfo.* files.
    # For each sfrinfo.* file, guided loop over particle ID.
    #   - a(sfrinfo) < a(fgal) is a sure thing.
    # Both files should be ranked with PID.
    #   - For a fast searching, should bin PID into bins
    # Each time sfrinfo.* finds a likely match for a PID, find
    #   - tmax(sfrinfo), tmax(fgal)
    # If tmax match, then YES!
    fformat={'names': ('atime', 'ID', 'alast', 'tmax', 'mstar'),
             'formats': ('f8', 'i8', 'f8', 'f8', 'f8')}
    # fsfrinfo = sfrinfobase + "sfrinfo.56"
    gal = sort(gal, order='ID')
    nparts = len(gal)
    aacc, alast, tmaxdiff, mstar = [-1]*nparts, [-1]*nparts, [10.]*nparts, [-1]*nparts
    for fi in range(NCPU):
        fsfrinfo = sfrinfobase + "sfrinfo." + str(fi)
        print "Reading: ", fsfrinfo
        acc = genfromtxt(fsfrinfo, dtype=fformat, usecols=(0,1,2,3,5))
        acc = sort(acc, order='ID')
        igal, nextPID = 0, gal['ID'][0]
        for iacc in range(len(acc)):
            if(acc['ID'][iacc] < nextPID): continue
            # Now acc.ID >= nextPID
            # while(gal['ID'][igal] == nextPID): # Start counting in the gal array
            while(acc['ID'][iacc] >= nextPID): # Start counting in the gal array  
                if(acc['ID'][iacc] == nextPID): # For each gal particle, see if a match is found
                    if(match_check(acc['tmax'][iacc], gal['Tmax'][igal], tol=tmaxdiff[igal])): # tmaxs match
                        if(acc['atime'][iacc] < gal['Age'][igal]+ATOL): # SF after accretion
                            aacc[igal] = acc['atime'][iacc]
                            alast[igal] = acc['alast'][iacc]
                            mstar[igal] = acc['mstar'][iacc]         
                            tmaxdiff[igal] = abs(acc['tmax'][iacc]-gal['Tmax'][igal])
                    # End (acc['ID'][iacc] == nextPID)
                igal += 1
                if(igal > nparts - 1): break # No particle to match!            
                # while(acc['ID'][iacc] >= nextPID):
                nextPID = gal['ID'][igal] # now updated.
    # Now we got aacc, alast, tmaxdiff.
    # Now Write:
    fout = open(outname, "w")
    fout.write("#PID a_form a_acc a_last Mass Tmax Mstar\n")
    for igal in range(len(gal)):
        if(alast[igal] < 0): gal['Tmax'][igal] *= -1
        line = "%8d %7.5f %7.5f % 7.5f %7.5e % 5.3f %6.3f\n" % \
        (gal['ID'][igal], gal['Age'][igal], aacc[igal], alast[igal], \
         gal['Mass'][igal], gal['Tmax'][igal], mstar[igal])
        fout.write(line)
    fout.close()
    return aacc, alast, tmaxdiff

def galsfh(galacc, fraction=False, lstyle="-"):
    galacc= sort(galacc, order='a_form')
    asf, mtot, mhot, mcold, mwhot, mwcold, mshot, mscold = \
        [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0]
    for igal in range(len(galacc)):
        asf.append(galacc['a_form'][igal])
        mhot.append(mhot[-1])
        mcold.append(mcold[-1])
        mwhot.append(mwhot[-1])
        mwcold.append(mwcold[-1])
        mshot.append(mshot[-1])
        mscold.append(mscold[-1])                
        # Tmax < 0 if last time it is a wind.
        if(galacc['a_last'][igal] > 0.): # Neither wind nor primodial accretion
            if(galacc['Tmax'][igal] > 5.5): mshot[-1] += galacc['Mass'][igal]
            else: mscold[-1] += galacc['Mass'][igal]            
        else:
            if(galacc['Tmax'][igal] > 5.5): mhot[-1] += galacc['Mass'][igal]
            elif(galacc['Tmax'][igal] > 0.0): mcold[-1] += galacc['Mass'][igal]
            elif(galacc['Tmax'][igal] > -5.5): mwcold[-1] += galacc['Mass'][igal]
            else: mwhot[-1] += galacc['Mass'][igal]
        mtot.append(mtot[-1]+galacc['Mass'][igal])
    if(fraction == False):
        p, = plt.plot(asf, mtot, linestyle=lstyle, color="black")
        plt.plot(asf, mcold, linestyle=lstyle, color="blue")
        plt.plot(asf, mhot, linestyle=lstyle, color="red")
        plt.plot(asf, mwcold, linestyle=lstyle, color="darkgreen")
        plt.plot(asf, mwhot, linestyle=lstyle, color="limegreen")
        plt.plot(asf, mshot, linestyle=lstyle, color="steelblue")
        plt.plot(asf, mscold, linestyle=lstyle, color="plum")
        l1 = plt.legend(["Total", "ColdP", "HotP", "ColdW", "HotW", "ColdS", "HotS"], fontsize=12)
        plt.gca().add_artist(l1)
        return p
    else:
        mcold = array(mcold)
        mhot = array(mhot)
        mwcold = array(mwcold)
        mwhot = array(mwhot)
        mshot = array(mshot)
        mscold = array(mscold)        
        mtot = array(mtot)
        p, = plt.plot(asf, (mtot/mtot[-1]), linestyle=lstyle, color="black")
        plt.plot(asf, (mcold/mtot), linestyle=lstyle, color="blue")
        plt.plot(asf, (mhot/mtot), linestyle=lstyle, color="red")
        plt.plot(asf, (mwcold/mtot), linestyle=lstyle, color="darkgreen")
        plt.plot(asf, (mwhot/mtot), linestyle=lstyle, color="limegreen")
        plt.plot(asf, (mscold/mtot), linestyle=lstyle, color="steelblue")
        plt.plot(asf, (mshot/mtot), linestyle=lstyle, color="plum")
        l1 = plt.legend(["Total", "ColdP", "HotP", "ColdW", "HotW", "ColdS", "HotS"], fontsize=12)
        plt.gca().add_artist(l1)
        return p

def merging_history(fgallst):
    # Show the merging history of a galaxy.
    # Basically plotting log(M*) - a_form 
    clrs = ["black", "blue", "red", "green", "orange"]
    ps = []
    flist = open(fgallst, "r")
    galacc= sort(galacc, order='a_form')
    for fi, line in enumerate(flist):
        fgalname = line.split()[0]
        print "Loading File: ", fgalname
        galacc = read_galacc_data(fgalname)
        # for igal in range(len(galacc)):
        #     if(galacc['a_last'][igal] < 0.):
        #         asf.append(galacc['a_form'][igal])
        #         tdiff = tcosmic(galacc['a_form'][igal]) - tcosmic(-galacc['a_last'][igal])
        #         trec.append(tdiff / 1.e9)
        # p, = plt.plot(asf, trec, ".", color = clrs[i], alpha=0.2, markersize=2)
        ps.append(p)
    flist.close()

def last_accretion(fgallst): # When do winds last time accreted?
    fig = plt.figure(1, figsize = (6, 6))
    flist = open(fgallst, "r")
    ls = ["-", "--"]
    clrs = ["blue", "red"]
    i = 0
    ps = []
    for line in flist:
        # gid = int(line.split()[0])
        # fgalname = "./p50n288zw/galacc_%d.stars" % (gid)
        fgalname = line.split()[0]
        print "Loading File: ", fgalname
        galacc = read_galacc_data(fgalname)
        asf, trec = [], []
        for igal in range(len(galacc)):
            if(galacc['a_last'][igal] < 0.):
                asf.append(galacc['a_form'][igal])
                tdiff = tcosmic(galacc['a_form'][igal]) - tcosmic(-galacc['a_last'][igal])
                trec.append(tdiff / 1.e9)
        p, = plt.plot(asf, trec, ".", color = clrs[i], alpha=0.2, markersize=2)
        ps.append(p)
        i = 1
    flist.close()
    plt.xlabel(r"$a_{form}$")
    plt.ylabel(r"$t_{rec} [Gyr]$")
    # plt.yscale("log")
    # plt.axis([0.0, 1.0, 1.e10, 5.e12])
    plt.legend(ps, ["Fiducial", "sw"], loc=4, fontsize=12)
    plt.show()

# gal = read_gal_data("gal_6662.stars")
# aacc, alast, tmaxdiff = load_sfrinfo(gal, fbase, "accgal_6662.stars")

def compile_galstar_info(fname):
    fgallst = open(fname, "r")
    for line in fgallst:
        gid = int(line.split()[0])
        print "Loading GID = ", gid
        fgalname = galbase+"gal_%d.stars" % (gid)
        foutname = galbase+"galacc_%d.stars" % (gid)        
        gal = read_gal_data(fgalname)
        aacc, alast, tmaxdiff = load_sfrinfo(gal, fbase, foutname)
    fgallst.close()

def show_starformation_history(fgallst, show_fraction=True):
    fig = plt.figure(1, figsize = (6, 6))
    flist = open(fgallst, "r")
    ls = ["-", "--"]
    i = 0
    ps = []
    for line in flist:
        # gid = int(line.split()[0])
        # fgalname = "./p50n288zw/galacc_%d.stars" % (gid)
        fgalname = line.split()[0]
        print "Loading File: ", fgalname
        galacc = read_galacc_data(fgalname)
        p = galsfh(galacc, fraction=show_fraction, lstyle=ls[i])
        ps.append(p)
        i = 1
    flist.close()
    plt.xlabel(r"$a(t)$")
    plt.ylabel(r"$M_*(a)$")
    # plt.yscale("log")
    plt.axis([0.0, 1.0, 1.e10, 5.e12])
    plt.legend(ps, ["Fiducial", "sw"], loc=4, fontsize=12)
