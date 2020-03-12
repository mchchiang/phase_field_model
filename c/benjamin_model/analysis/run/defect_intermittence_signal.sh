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
    echo "usage: defect_intermittence_signal.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=36 # 100
tstart=1000000
tend=21000000
tinc=1000
defect_thres=0.1 #0.05
jump_thres=2
time_thres=100000

intermit_exe="../bin/exe/defect_intermittence_signal"

max_jobs=8 # 8
cmd=()
jobid=0

while (( $(bc <<< "$d < $d_end") ))
do
    #in_path="${in_dir}/d_${d}/neighdiff"
    in_path="${in_dir}/d_${d}/neigh_delaunay"
    if [ -d $in_path ]; then
	out_path="${out_dir}/d_${d}/intermit/"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
		neighdiff_file="${in_path}/neighdiff_${name}.dat"
		if [ -f $neighdiff_file ]; then
		    out_file="${out_path}/intermit_${name}.dat"
		    echo "Doing d = $d Pe = $pe run = $run"
		    cmd[$jobid]="$intermit_exe $tstart $tend $tinc $defect_thres $time_thres $neighdiff_file $out_file"
		    jobid=$(bc <<< "$jobid + 1")
		fi
	    done
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done

# Parallel runs

total_jobs=$jobid
jobid=0

while (( $(bc <<< "$jobid < $total_jobs") ))
do
    for (( i=0; i<$max_jobs && $jobid < $total_jobs; i++))
    do
	echo "${cmd[jobid]} &"
	${cmd[jobid]} &
	jobid=$(bc <<< "$jobid + 1")
    done
    wait
done
