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
    echo "usage: msd_tavg.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

msd_exe="../bin/exe/msd"

N=100
tstart=0
tend=20000000
tinc=1000

max_jobs=10
nthreads=1
cmd=()
jobid=0
export OMP_NUM_THREADS=$nthreads

while (( $(bc <<< "$d < $d_end") ))
do
    pe=$(python -c "print '%.3f' % ($pe_start)")
    while (( $(bc <<< "$pe < $pe_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
	    pos_file="${in_dir}/position/pos_${name}.dat"
	    params_file="${in_dir}/siminfo/params_${name}.txt"
	    lx=$(grep 'lx = ' $params_file | awk '{print $3}')
	    ly=$(grep 'ly = ' $params_file | awk '{print $3}')
	    msd_file="${out_dir}/msd_${name}.dat"
	    cmd[$jobid]="$msd_exe $N $lx $ly $tstart $tend $tinc $pos_file $msd_file"
	    jobid=$(bc <<< "$jobid + 1")
	done
	pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
    done
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
