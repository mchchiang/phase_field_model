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

if [ "$#" != 10 ]; then
    echo "usage: check_topocharge.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir"
    exit 1
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

topo_py="../src/check_topocharge.py"

N=100 #100
tstart=0
tend=21000000
tinc=1000

max_jobs=8 # 8
cmd=()
jobid=0


while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/neigh_delaunay/"
    if [ -d $in_path ]; then
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
		neigh_file="${in_path}/neigh_${name}.dat"
		if [ -f $neigh_file ]; then
		    cmd[$jobid]="python $topo_py $N $neigh_file"
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
