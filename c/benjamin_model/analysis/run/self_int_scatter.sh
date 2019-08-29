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
    echo "usage: int_scatter_fn.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

sisf_exe="../bin/exe/self_int_scatter"

N=100
tstart=0
tend=20000000
tinc=1000
tshiftend=20000000 # 20000000
r0=8
nqvec=10

max_jobs=1
nthreads=4
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
	    pos_bulk_file="${in_dir}/position/pos-bulk_${name}.dat"
	    params_file="${in_dir}/siminfo/params_${name}.txt"
	    lx=$(grep 'lx = ' $params_file | awk '{print $3}')
	    ly=$(grep 'ly = ' $params_file | awk '{print $3}')
	    sisf_file="${out_dir}/sisf_${name}_r0_${r0}_ts_${tshiftend}.dat"
	    cmd[$jobid]="$sisf_exe $N $lx $ly $nqvec $r0 $tstart $tend $tinc $tshiftend $pos_file $pos_bulk_file $sisf_file"
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
