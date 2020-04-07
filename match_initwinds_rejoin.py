from numpy import genfromtxt, log10, sqrt
from cosmology import tcosmic
from astroconst import pc, ac
import phew
import matplotlib.pyplot as plt
import sys

# Plan:
# Select winds first, match later.
# VERSION 1.0:
# Proto-type. Need initwinds.sorted and rejoin.sorted first

# SEARCH "------" FOR MAIN()

# In z?/:
# tabinit["atime","Mstar","PhEWKey","T_a","Vinit"]
# tabrejoin["atime","Vinit","PhEWKey","M_c","vrel", "cs_a", "rho_a","flag"]
# snapshot_???.phews (idx, PhEWKey, HID) >> Mvir, Msub, Rvir
#   - $SCI-PHEW/select/find_host_haloes_for_phews.py
# Out:
# For all winds launched within a redshift window:
# a_init, Mvir, Msub, Rvir, Vinit, T_a (<- init:rejoin ->), a_rejoin, M_c, Mach, PhEWKey

__all__ = ['load_tables', 'match_tables']

errormsg = "Usage: "+sys.argv[0]+" ncpu redshift(0,1,2) modelname"
if(len(sys.argv) != 4):
    print errormsg
    sys.exit(1)
else:
    MODEL = sys.argv[1]
    zstr = sys.argv[2]
    snapstr = sys.argv[3]

# MODEL = "l25n144-phew"
# zstr = "z0"
# snapstr = "100"

def load_tables():
    print "Loading tables ... "
    FieldsInit = phew.PhEWFields("fields_initwinds.dat")
    filename = "/proj/shuiyao/"+MODEL+"/WINDS/"+zstr+"/initwinds.sorted"
    info_fields = FieldsInit.get_field_info(["atime","Mstar","PhEWKey","T_a","Vinit"], verbose=True)
    tabinit = phew.read_phew_fields(filename, info_fields)

    FieldsRejoin = phew.PhEWFields("fields_rejoin.dat")
    filename = "/proj/shuiyao/"+MODEL+"/WINDS/"+zstr+"/rejoin.sorted"
    info_fields = FieldsRejoin.get_field_info(["atime","Vinit","PhEWKey","M_c","vrel", "cs_a", "rho_a","flag"], verbose=True)
    tabrejoin = phew.read_phew_fields(filename, info_fields)

    filename = "/proj/shuiyao/"+MODEL+"/snapshot_"+snapstr+".phews"
    tabhalo = genfromtxt(filename, names=True)
    return tabinit, tabrejoin, tabhalo

def match_tables(tabi, tabr, tabh, filename):
    key_to_idx_r = dict()
    key_to_idx_h = dict()
    for i in range(len(tabr)):
        key_to_idx_r[tabr['PhEWKey'][i]] = i
    for i in range(len(tabh)):
        key_to_idx_h[tabh['PhEWKey'][i]] = i
    print "Matching Tables: ..."
    print "Launched: ", len(tabi)
    print "Halo found: ", len(tabh)
    print "Rejoined: ", len(tabr)
    counth, countr = 0, 0
    fout = open(filename, "w")
    fout.write("#a_i Vinit T_a a_rejoin M_c Mach LogMvir LogMsub Rvir PhEWKey\n")
    for i in range(len(tabi)):
        key = tabi['PhEWKey'][i]
        if(key in key_to_idx_r): countr += 1
        if(key in key_to_idx_h): counth += 1        
        # Sometimes there are very few PhEW that rejoined.
        # if(key in key_to_idx_r and key in key_to_idx_h):
        if(key in key_to_idx_h):
            fout.write("%7.5f %6.1f %5.3f " %
                       (tabi[i]['atime'], tabi[i]['Vinit'], log10(tabi[i]['T_a'])))
            if(key in key_to_idx_r):
                rejoin = tabr[key_to_idx_r[key]]
                fout.write("%7.5f %5.3f %5.2f " %
                           (rejoin['atime'], rejoin['M_c'], rejoin['vrel']/rejoin['cs_a']))
            else:
                fout.write("%7.5f %5.3f %5.2f " % (-1.0, -1.0, -1.0))
            halo = tabh[key_to_idx_h[key]]
            fout.write("%6.3f %6.3f %6.1f " % (halo['LogMvir'], halo['LogMsub'], halo['Rvir']))
            fout.write("%8d\n" % (key))
    fout.close()
    print "Launched - Rejoin pairs: ", countr
    print "Launched - Halo pairs: ", counth

def load_matched_data(filename):
    return genfromtxt(filename, dtype='i8')

def show_matched_data(tab_init, tab_rejoin, idxmap, nskip=1):
    iskip = 0
    flag_count = [0, 0, 0, 0, 0, 0]
    xarr, yarr = [], []
    for i in range(tab_init.size):
        if(iskip >= nskip):
            idx = idxmap[i]
            if(idx != -1):
                # print tabinit['atime'][i], tabrejoin['atime'][idx], tabinit['T_a'][i], tab_rejoin['M_c'][idx]
                # xarr.append(log10(tabinit['T_a'][i]))
                xarr.append(tabinit['Mstar'][i])
                # yarr.append((tcosmic(tabrejoin['atime'][idx])-tcosmic(tabinit['atime'][i]))/1.e6) # Myr
                yarr.append(tabrejoin['M_c'][idx])
                # if(yarr[-1] > 1000.):
                #     print tabinit['ID'][i], tabinit['Vinit'][i], tabrejoin['Vinit'][idx]
                if(yarr[-1] < 1.0):
#                    print "Mach = %5.3f; T_a = %3.1f; Vi = %5.1f; flag = %d" % (tabinit[i]['Vinit'] * 1.e5 / sqrt(pc.k * tabinit[i]['T_a'] / (0.6 * pc.mh)), log10(tabinit[i]['T_a']), tabinit[i]['Vinit'], tabrejoin[idx]['flag'])
                    flag_count[tabrejoin['flag'][idx]] += 1
                iskip = 0
        iskip += 1
    print flag_count
    plt.plot(xarr, yarr, "b.", alpha=0.2, markersize=3)
    plt.show()

# ------
foutname = "/scratch/shuiyao/scidata/newwind/"+MODEL+"/phewsinfo."+zstr
tabinit, tabrejoin, tabhalo = load_tables()
match_tables(tabinit, tabrejoin, tabhalo, foutname)
# rejoin_idx = load_matched_data("winds.dat")
# show_matched_data(tabinit, tabrejoin, rejoin_idx)

# Rejoin Flags:
# 0/1: Mach number small enough .OR. Velocity difference small enough
#  - 0: Temperature difference LARGE
#  - 1: Temperature difference SMALL


# OBSOLETE: Now we use hash map to match tables
def match_tables_old():    
    # Use Vinit and ID to match
    N_init, N_rejoin = len(tabinit), len(tabrejoin)
    i_init, i_rejoin = 0, 0
    fout = open("winds.dat", "w")
    buf_rejoin = []
    i_rejoin_start = 0
    thisIDrejoin = tabrejoin['PhEWKey'][i_rejoin]
    while(i_init < N_init):
        # 1. ID < ID
        if(tabinit['PhEWKey'][i_init] < thisIDrejoin): # Not rejoined
            i_init += 1
            fout.write("-1\n")
            continue
        # 2. ID > ID (end)
        if(tabinit['PhEWKey'][i_init] > thisIDrejoin and i_rejoin == N_rejoin):
            i_init += 1
            fout.write("-1\n")
            continue
        # ID >= ID;
        # 1. ID > ID: 
        if(tabinit['PhEWKey'][i_init] > thisIDrejoin): # Fetch for the new ID
            while(i_rejoin < N_rejoin and tabinit['PhEWKey'][i_init] > tabrejoin['PhEWKey'][i_rejoin]):
                i_rejoin += 1
            # now ID <= ID
            if(i_rejoin == N_rejoin): continue # All the rest hasn't rejoined
            buf_rejoin = []
            i_rejoin_start = i_rejoin
            thisIDrejoin = tabrejoin['PhEWKey'][i_rejoin]
            # Store all rejoin info for this new ID >>>
            while(i_rejoin < N_rejoin and tabrejoin['PhEWKey'][i_rejoin] == thisIDrejoin):
                buf_rejoin.append(tabrejoin[i_rejoin])
                i_rejoin += 1
            # <<<
        else: # ID == ID
            match = False
            for i in range(len(buf_rejoin)):
                if(tabinit['Vinit'][i_init] == buf_rejoin[i]['Vinit']):
                    if(tabinit['atime'][i_init] <= buf_rejoin[i]['atime']):
                        fout.write("%d\n" % (i_rejoin_start + i))
                        # print tabinit['atime'][i_init], tabrejoin['atime'][i_rejoin_start+i], tabinit['PhEWKey'][i_init], tabrejoin['PhEWKey'][i_rejoin_start+i]
                        match = True
                        break
            if(match == False): fout.write("-1\n")
            i_init += 1
    fout.close()

