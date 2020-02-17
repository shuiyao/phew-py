from numpy import genfromtxt, log10, sqrt
from cosmology import tcosmic
from astroconst import pc, ac
import phew
import matplotlib.pyplot as plt

__all__ = ['load_tables', 'match_tables']

def load_tables():
    print "Loading tables ... "
    FieldsInit = phew.PhEWFields("fields_initwinds.dat")
    filename = "/proj/shuiyao/m6n64beta5/WINDS/initwinds.sorted"
    info_fields = FieldsInit.get_field_info(["atime","ID","T_a","Vinit"], verbose=True)
    tabinit = phew.read_phew_fields(filename, info_fields)

    FieldsRejoin = phew.PhEWFields("fields_rejoin.dat")
    filename = "/proj/shuiyao/m6n64beta5/WINDS/rejoin.sorted"
    info_fields = FieldsRejoin.get_field_info(["atime","Vinit","ID","M_c","vrel", "cs_a", "rho_a","flag"], verbose=True)
    tabrejoin = phew.read_phew_fields(filename, info_fields)
    return tabinit, tabrejoin

def match_tables():    
    # Use Vinit and ID to match
    N_init, N_rejoin = len(tabinit), len(tabrejoin)
    i_init, i_rejoin = 0, 0
    fout = open("winds.dat", "w")
    buf_rejoin = []
    i_rejoin_start = 0
    thisIDrejoin = tabrejoin['ID'][i_rejoin]
    while(i_init < N_init):
        # 1. ID < ID
        if(tabinit['ID'][i_init] < thisIDrejoin): # Not rejoined
            i_init += 1
            fout.write("-1\n")
            continue
        # 2. ID > ID (end)
        if(tabinit['ID'][i_init] > thisIDrejoin and i_rejoin == N_rejoin):
            i_init += 1
            fout.write("-1\n")
            continue
        # ID >= ID;
        # 1. ID > ID: 
        if(tabinit['ID'][i_init] > thisIDrejoin): # Fetch for the new ID
            while(i_rejoin < N_rejoin and tabinit['ID'][i_init] > tabrejoin['ID'][i_rejoin]):
                i_rejoin += 1
            # now ID <= ID
            if(i_rejoin == N_rejoin): continue # All the rest hasn't rejoined
            buf_rejoin = []
            i_rejoin_start = i_rejoin
            thisIDrejoin = tabrejoin['ID'][i_rejoin]
            # Store all rejoin info for this new ID >>>
            while(i_rejoin < N_rejoin and tabrejoin['ID'][i_rejoin] == thisIDrejoin):
                buf_rejoin.append(tabrejoin[i_rejoin])
                i_rejoin += 1
            # <<<
        else: # ID == ID
            match = False
            for i in range(len(buf_rejoin)):
                if(tabinit['Vinit'][i_init] == buf_rejoin[i]['Vinit']):
                    if(tabinit['atime'][i_init] <= buf_rejoin[i]['atime']):
                        fout.write("%d\n" % (i_rejoin_start + i))
                        # print tabinit['atime'][i_init], tabrejoin['atime'][i_rejoin_start+i], tabinit['ID'][i_init], tabrejoin['ID'][i_rejoin_start+i]
                        match = True
                        break
            if(match == False): fout.write("-1\n")
            i_init += 1
    fout.close()

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
                xarr.append(log10(tabinit['T_a'][i]))
                yarr.append((tcosmic(tabrejoin['atime'][idx])-tcosmic(tabinit['atime'][i]))/1.e6) # Myr
                # if(yarr[-1] > 1000.):
                #     print tabinit['ID'][i], tabinit['Vinit'][i], tabrejoin['Vinit'][idx]
                if(yarr[-1] < 1.0):
                    # print tabinit[i], tabrejoin[idx]
#                    print "Mach = %5.3f; T_a = %3.1f; Vi = %5.1f; flag = %d" % (tabinit[i]['Vinit'] * 1.e5 / sqrt(pc.k * tabinit[i]['T_a'] / (0.6 * pc.mh)), log10(tabinit[i]['T_a']), tabinit[i]['Vinit'], tabrejoin[idx]['flag'])
                    flag_count[tabrejoin['flag'][idx]] += 1
                iskip = 0
        iskip += 1
    print flag_count
    plt.plot(xarr, yarr, "b.", alpha=0.2, markersize=3)
    plt.show()

# tabinit, tabrejoin = load_tables()
match_tables()
rejoin_idx = load_matched_data("winds.dat")
show_matched_data(tabinit, tabrejoin, rejoin_idx)

# Rejoin Flags:
# 0/1: Mach number small enough .OR. Velocity difference small enough
#  - 0: Temperature difference LARGE
#  - 1: Temperature difference SMALL
