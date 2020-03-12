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
    echo "usage: compress_field.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=400 #100
tstart=0
tend=21000000
tinc=1000

max_jobs=8 # 8
cmd=()
jobid=0
export XZ_OPT=-9
while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/field/"
    if [ -d $in_path ]; then
	out_path="${out_dir}/d_${d}/field_xz/"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
		field_file="${in_path}/field_${name}.dat.0"
		if [ -f $field_file ]; then 
		    field_xz="${out_path}/field_${name}.tar.xz"
		    cmd[jobid]="files=\$(ls ${in_path}/field_${name}.dat.* | xargs -n 1 basename) && tar -Jcvf ${field_xz} -C ${in_path} \$files"
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
	eval "${cmd[jobid]} &"
	jobid=$(bc <<< "$jobid + 1")
    done
    wait
done
