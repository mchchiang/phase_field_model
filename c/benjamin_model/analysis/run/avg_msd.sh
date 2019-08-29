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
    echo "usage: msd_tavg.sh d_start d_end d_inc pe_start pe_end pe_inc in_dir out_dir"
    exit 1
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

avg_py="../src/AverageMultiFiles.py"

N=100
tstart=0
tend=20000000
tinc=1000


while (( $(bc <<< "$d < $d_end") ))
do
    pe=$(python -c "print '%.3f' % ($pe_start)")
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	name="cell_N_${N}_d_${d}_Pe_${pe}"
	msd="${in_dir}/msd_${name}"
	avg_msd_file="${in_dir}/msd_${name}_avg.dat"
	python $avg_py 0 1 -1 -1 $avg_msd_file "${msd}_run_"*.dat
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
