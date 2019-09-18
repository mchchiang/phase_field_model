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
    echo "usage: prob_hexatic.sh d_start d_end d_inc pe_start pe_end pe_inc in_dir out_dir"
    exit 1
fi

prob_py="../src/Distribution.py"

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=100 #100
run=2

tstart=10000000
tend=20000000
tinc=1000

min=0
max=1
bin_size=0.02

while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/hexatic/"
    pe=$(python -c "print '%.3f' % ($pe_start)")
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
	hexatic_file="${in_path}/hexatic_${name}.dat"
	if [ -f $hexatic_file ]; then
	    out_path="${out_dir}/d_${d}/hexatic/prob/t_${tstart}-${tend}/"
	    if [ ! -d $out_path ]; then
		mkdir -p $out_path
	    fi

	    # Filter for selected time frames
	    data_file="${out_path}/hexatic-filtered_${name}.dat"
	    awk -v ts=${tstart} -v te=${tend} '{if($1>=ts&&$1<=te){print}}' $hexatic_file > $data_file
	    prob_file="${out_path}/prob-hexatic-mag-avg_${name}.dat"
	    python $prob_py 3 $min $max $bin_size $data_file $prob_file
	    rm $data_file
	fi
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
