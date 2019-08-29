#!/bin/bash
d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
in_dir=$7
out_dir=$8

if [ "$#" != 8 ]; then
    echo "usage: msd_exp.sh d_start d_end d_inc pe_start pe_end pe_inc in_dir out_dir"
    exit 1
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=100
tstart=0
tend=20000000
tinc=1000
tshiftend=20000000 # 20000000

msd_exp_file=${out_dir}/msd-exp_cell_N_${N}_d_${d_start}-${d_end}_Pe_${pe_start}-${pe_end}.dat
> $msd_exp_file

a=1.0
b=1.0
xmin=2.0
xmax=3.0
fitgp="msd_fit.gp"
fitlog="msd_fit.log"

while (( $(bc <<< "$d < $d_end") ))
do
    pe=$(python -c "print '%.3f' % ($pe_start)")
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	name="cell_N_${N}_d_${d}_Pe_${pe}"
	msd_file="${in_dir}/d_${d}/msd/msd-tavg_${name}_ts_${tshiftend}_avg.dat"
	if [ -f $msd_file ]; then
	    > $fitgp # Clear the fit file and fit log
	    > $fitlog
	    echo "f(x) = a*x+b
a=${a}
b=${b}
set fit logfile '${fitlog}'
fit [${xmin}:${xmax}] f(x) '${msd_file}' u (log10(\$1*0.00005)):(log10(\$2/12**2)) via a,b
" > $fitgp
	    echo "Fitting MSD exponent for d = ${d} Pe = ${pe}"
	    gnuplot $fitgp
	    data_a=$(grep -E "a +[=] [-0-9.]* +[+][/][-]" $fitlog | awk '{print $3,$5}')
	    data_b=$(grep -E "b +[=] [-0-9.]* +[+][/][-]" $fitlog | awk '{print $3,$5}')
	    echo $d $pe $data_a $data_b >> $msd_exp_file
	fi
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done

