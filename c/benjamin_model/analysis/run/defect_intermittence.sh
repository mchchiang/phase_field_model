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
    echo "usage: defect_intermittence.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=400 # 100
#run=1
tstart=1000000
tend=21000000
tinc=1000
defect_thres=0.05 #0.05
jump_thres=2
time_thres=100000

intermit_exe="../bin/exe/defect_intermittence"

for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
do
    intermit_file="${out_dir}/intermit_cell_N_${N}_run_${run}_t_${time_thres}.dat"
    #intermit_file="intermit_cell_N_${N}_run_${run}_t_${time_thres}.dat"
    > $intermit_file
    d=$(python -c "print '%.3f' % ($d_start)")
    while (( $(bc <<< "$d < $d_end") ))
    do
	#in_path="${in_dir}/d_${d}/neighdiff"
	in_path="${in_dir}/d_${d}/neigh_delaunay"
	if [ -d $in_path ]; then
	    #out_path="${out_dir}/d_${d}/neighdiff"
	    #if [ ! -d $out_path ]; then
	    #    mkdir -p $out_path
	    #fi
	    pe=$(python -c "print '%.3f' % ($pe_start)")
	    while (( $(bc <<< "$pe < $pe_end") ))
	    do
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
		neighdiff_file="${in_path}/neighdiff_${name}.dat"
		if [ -f $neighdiff_file ]; then
		    echo "Doing d = $d Pe = $pe run = $run"
		    data=$($intermit_exe $tstart $tend $tinc $defect_thres $jump_thres $time_thres $neighdiff_file)
		    echo "$d $pe $data" >> $intermit_file
		fi
		pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	    done
	    echo >> $intermit_file
	fi
	d=$(python -c "print '%.3f' % ($d + $d_inc)")
    done
done
