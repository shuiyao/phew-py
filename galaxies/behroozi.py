import ioformat
import matplotlib.pyplot as plt
from scipy import array
from scipy.interpolate import interp1d

# c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat
# c_smmr_then_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat

# /scratch/shuiyao/sci/evolve/ghistory/data_behroozi13/sm/sm_hist_rel_11.0.dat

zstrs = ["z0.10", "z1.00", "z2.00", "z3.00", "z4.00"]
smmr_name_prefix = "c_smmr_"
smmr_name_suffix = "_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat"
SMMRBASE = "/scratch/shuiyao/sci/REFERENCES/behroozi13/smmr/"
SMBASE = "/scratch/shuiyao/sci/REFERENCES/behroozi13/sm/"

def read_smmr(z, mhthen=False):
    fname = SMMRBASE
    fname += smmr_name_prefix
    if(mhthen == True): fname += "then_"
    i = (int)(z)
    fname += zstrs[i] + smmr_name_suffix
    print "Reading: ", fname
    mh, dm, e1, d2 = ioformat.rcol(fname, [0,1,2,3], linestart=1)
    return mh, dm, e1, d2

mhstrs = ["11.0", "12.0", "13.0", "14.0", "15.0"]
sm_name_prefix = "sm_hist_rel_"
sm_name_suffix = ".dat"

def read_sm(mh):
    fname = SMBASE
    fname += sm_name_prefix
    i = (int)(mh - 11.0)
    fname += mhstrs[i] + sm_name_suffix
    print "Reading: ", fname
    a, frac, e1, e2 = ioformat.rcol(fname, [0,1,2,3])
    return a, frac, e1, e2

def plot_smmr(z):
    mh, dm, e1, e2 = read_smmr(z)
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.errorbar(mh, dm, yerr=[e1, e2], color="black", fmt='o')
    plt.show()

def ms2mh(z, kind='linear'):
    mh, dm, e1, e2 = read_smmr(z)
    ms=array(mh)+array(dm)
    f = interp1d(ms, mh, kind=kind, fill_value='extrapolate')
    return f
