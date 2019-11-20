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
    echo "usage: avg_msd_tavg.sh d_start d_end d_inc pe_start pe_end pe_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

avg_py="../src/AverageMultiFiles.py"

N=100

while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/msd/t_1000000-21000000/"
    if [ -d $in_path ]; then
	out_path="${out_dir}/d_${d}/msd/t_1000000-21000000/"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    name="cell_N_${N}_d_${d}_Pe_${pe}"
	    msd="${in_path}/msd-tavg_${name}"
	    if [ -f ${msd}_run_1.dat ]; then
		echo "Doing N = ${N} d = ${d} Pe = ${pe}"
		avg_msd_file="${out_path}/msd-tavg_${name}_avg.dat"
		python $avg_py 0 1 -1 -1 $avg_msd_file "${msd}_run_"*.dat
		awk '{print int($1),$2,$3,$4}' $avg_msd_file > ${avg_msd_file}.tmp
		mv ${avg_msd_file}.tmp ${avg_msd_file}
	    fi
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
