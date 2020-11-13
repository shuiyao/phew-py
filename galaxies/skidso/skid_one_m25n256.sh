#!/usr/bin/bash

#Usage: bash skid_one.sh snapnum redshift runname

database="/proj/shuiyao/"
runname=$3
prefix=$database$runname/"snapshot_" # for binary naming scheme prefiXz$zSuffix
reducedir=$database$runname/

# NOTE: DON'T FORGET TO CHANGE Z
if [ ! $2 ] 
then
    echo "Not enough paramters."
    exit 0
fi
num=$1
if [ ${#num} -eq 1 ]
then
    num="00"$num
elif [ ${#num} -eq 2 ]
then
    num="0"$num
elif [ ${#num} -eq 3 ]
then
    num=$num
else
    echo "Wrong Parameter."
    exit 0
fi
z=$2

suffix=".bin"      # e.g. cosmoL144z3.0.bin

binfile=$prefix$num$suffix
auxfile=$prefix$num".aux"
outputfile=$reducedir"/gal_z"$num
#mkdir $reducedir
echo "SKID <- Binary: "$binfile
echo "outputfile: "$outputfile

#simulation parameters - likely to change
omega0=0.30	# value of omega (choose one to give right high-z H_0)
Lambda=0.70	# value of cosmological constant
h0=0.68           # Hubble constant in 100 km/sec/Mpc units
omegabh2=0.048 # value of omegab (used to determine threshold density)
#eps=7.0e-5   # gravitational softening parameter in system units
eps=3.5e-5   # p50n576
nsm=32		# number of particles for smoothing
h=2.8944        # Hubble constant in simulation units
p=1		# period of simulation cube
nmin=8		# minimum number of members of SKID groups
omegab=`echo "scale=5;$omegabh2/$h0/$h0" | bc`
tau=`python -c "print 2.*$eps"`           # linking length
echo "eps, omegab = "$eps", "$omegab
dens=`python -c "print $omegab*1000."`
#fofll=0.14/$nside/2.**(1/3.) 

# OBSOLETE
# you are less likely to want to change these
# munit=2.77579e11*$lbox*$lbox*$lbox/$h0
# lmpc=$lbox/$h0
# nside=576     #cube of this is number of particles (in each of gas and dark)
# lbox=50.0         # size of box in h^-1 Mpc
# temp=30000	# temperature threshold
# npart=2*$nside*$nside*$nside

echo "skid -twophase -sfrfile "$auxfile" -tau "$tau" -e "$eps" -s "$nsm" -d "$dens" -H "$h" -O "$omega0" -Lambda "$Lambda" -z "$z" -p "$p" -m "$nmin" -stats -o "$outputfile" < "$binfile

/home/shuiyao/bin/skid -twophase -sfrfile $auxfile -tau $tau -e $eps -s $nsm -d $dens -H $h -O $omega0 -Lambda $Lambda -z $z -p $p -m $nmin -stats -o $outputfile < $binfile

# echo "skid -twophase -tau "$tau" -e "$eps" -s "$nsm" -d "$dens" -t "$temp" -H "$h" -O "$omega0" -Lambda "$Lambda" -z "$z" -p "$p" -m "$nmin" -I "$binfile" -stats -o "$outputfile

# skid -twophase -tau $tau -e $eps -s $nsm -d $dens -t $temp -H $h -O $omega0 -Lambda $Lambda -z $z -p $p -m $nmin -I $binfile -stats -o $outputfile
