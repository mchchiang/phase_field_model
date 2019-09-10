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
    echo "usage: avg_eccent.sh d_start d_end d_inc pe_start pe_end pe_inc in_dir out_dir"
    exit 1
fi

time_avg_py="../src/TimeAverage.py"

d=$(python -c "print '%.3f' % ($d_start)")
d_old=$d
pe=$(python -c "print '%.3f' % ($pe_start)")

N=36 #100
run=1

eccent_tstart=0
eccent_tend=20000000
tstart=10000000
tend=20000000
tinc=1000

eccent_avg_file="${out_dir}/eccent_cell_N_${N}_d_${d_start}-${d_end}_Pe_${pe_start}-${pe_end}_t_${tstart}-${tend}_run_${run}.dat"
> $eccent_avg_file

while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/eccent/"
    pe=$(python -c "print '%.3f' % ($pe_start)")
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
	eccent_file="${in_path}/eccent_${name}.dat"
#	eccent_file="${in_path}/eccent_${name}_t_${eccent_tstart}-${eccent_tend}.dat"
	if [ -f $eccent_file ]; then
	    if [ $d != $d_old ]; then
		echo "" >> $eccent_avg_file
		d_old=$d
	    fi
	    out_path="${out_dir}/d_${d}/eccent/"
	    if [ ! -d $out_path ]; then
		mkdir -p $out_path
	    fi
	    avg_file="${out_path}/eccent_${name}_t_${tstart}-${tend}_avg.dat"
	    python $time_avg_py 0 1 $tstart $tend $tinc $eccent_file $avg_file
	    data=$(cat $avg_file)
	    echo "$d $pe $data" >> $eccent_avg_file 
	    rm $avg_file
	fi
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
