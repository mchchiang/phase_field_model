#!/bin/bash
d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
run_start=$7
run_end=$8
run_inc=$9
qsub_index_file=${10}

if [ "$#" != 10 ]; then
    echo "usage: phase_field_multi_init.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc qsub_index_file"
    exit 1
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=100

> $qsub_index_file # Clear the index file

while (( $(bc <<< "$d < $d_end") ))
do
    pe=$(python -c "print '%.3f' % ($pe_start)")
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    dir="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
	    echo $dir >> $qsub_index_file
	done
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
