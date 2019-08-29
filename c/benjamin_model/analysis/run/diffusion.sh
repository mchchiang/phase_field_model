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
    echo "usage: diffusion.sh d_start d_end d_inc pe_start pe_end pe_inc in_dir out_dir"
    exit 1
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

msd_exe="../bin/exe/msd_time_avg.2"

N=100
tstart=0
tend=20000000
tinc=1000
tshiftend=20000000 # 20000000

deff_file=${out_dir}/deff_cell_N_${N}_d_${d_start}-${d_end}_Pe_${pe_start}-${pe_end}.dat
> $deff_file

while (( $(bc <<< "$d < $d_end") ))
do
    pe=$(python -c "print '%.3f' % ($pe_start)")
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	name="cell_N_${N}_d_${d}_Pe_${pe}"
	msd_file="${in_dir}/d_${d}/msd/msd-tavg_${name}_ts_${tshiftend}_avg.dat"
	echo $msd_file
	if [ -f $msd_file ]; then
	    in_data=$(tail -n1 $msd_file)
	    tstep=$(echo $in_data | awk '{print $1}')
	    msd_val=$(echo $in_data | awk '{print $2}')
	    deff=$(python -c "print '%.10g' % ($msd_val / (4.0*$tstep))")
	    echo $d $pe $deff >> $deff_file
	fi
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done

