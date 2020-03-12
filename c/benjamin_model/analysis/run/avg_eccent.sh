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
in_dir=${10}
out_dir=${11}

if [ "$#" != 11 ]; then
    echo "usage: avg_eccent.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

time_avg_py="../src/TimeAverage.py"
multi_avg_py="../src/AverageMultiFiles.py"

d=$(python -c "print '%.3f' % ($d_start)")
d_old=$d
pe=$(python -c "print '%.3f' % ($pe_start)")

N=36 #100

tstart=1000000
tend=21000000
tinc=1000

eccent_overall_file="${out_dir}/eccent_cell_N_${N}_d_${d_start}-${d_end}-${d_inc}_Pe_${pe_start}-${pe_end}-${pe_inc}_t_${tstart}-${tend}.dat"
> $eccent_overall_file

while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/eccent/"
    pe=$(python -c "print '%.3f' % ($pe_start)")
    
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	name="cell_N_${N}_d_${d}_Pe_${pe}"
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    eccent_file="${in_path}/eccent_${name}_run_${run}.dat"
	
	    if [ -f $eccent_file ]; then
		echo "Doing d = $d pe = $pe run = $run"
		out_path="${out_dir}/d_${d}/eccent/t_${tstart}-${tend}"
		if [ $d != $d_old ]; then
		    echo "" >> $eccent_overall_file
		    d_old=$d
		fi
		if [ ! -d $out_path ]; then
		    mkdir -p $out_path
		fi
		time_avg_file="${out_path}/eccent_${name}_run_${run}_avg.dat"
		python $time_avg_py 0 1 $tstart $tend $tinc $eccent_file $time_avg_file
	    fi
	done
	if [ -f "${out_path}/eccent_${name}_run_1_avg.dat" ]; then
	    # Average over differen runs
	    avg_file="${out_path}/eccent_${name}_avg.dat"
	    python $multi_avg_py -1 0 2 -1 $avg_file "${out_path}/eccent_${name}"_run_*_avg.dat
	    data=$(cat $avg_file)
	    echo "$d $pe $data" >> $eccent_overall_file 
	    rm "${out_path}/eccent_${name}"_run_*_avg.dat
	    rm $avg_file
	fi
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
