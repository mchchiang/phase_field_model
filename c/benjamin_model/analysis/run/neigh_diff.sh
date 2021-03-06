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
    echo "usage: neigh_diff.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=400 # 100
tstart=0
tend=21000000

max_jobs=8 # 8
cmd=()
jobid=0

neigh_exe="../bin/exe/neigh_diff"

while (( $(bc <<< "$d < $d_end") ))
do
    #in_path="${in_dir}/d_${d}/neighbour"
    in_path="${in_dir}/d_${d}/neigh_delaunay"
    #in_path=$in_dir
    if [ -d $in_path ]; then
	out_path="${out_dir}/d_${d}/neigh_delaunay"
	#out_path=$out_dir
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
		neigh_file="${in_path}/neigh_${name}.dat"
		if [ -f $neigh_file ]; then
		    echo "Doing d = $d Pe = $pe run = $run"
		    neighdiff_file="${out_path}/neighdiff_${name}.dat"
		    cmd[$jobid]="$neigh_exe $N $tstart $tend $neigh_file $neighdiff_file"
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

