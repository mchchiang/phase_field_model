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
    echo "usage: overlap_defect_corr.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

olap_exe="../bin/exe/overlap_defect_corr"

N=100 #100
tstart=1000000
tend=21000000
tinc=1000
ftinc=10000

max_jobs=8 # 8
cmd=()
jobid=0

while (( $(bc <<< "$d < $d_end") ))
do
    #in_path="${in_dir}/d_${d}/"
    in_path="${in_dir}"
    if [ -d $in_path ]; then
	#out_path="${out_dir}/d_${d}/local_overlap/"
	out_path="${out_dir}/local_overlap_2"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
		pos_file="${in_path}/position/pos_${name}.dat"
		if [ -f $pos_file ]; then
		    neigh_file="${in_path}/neigh_delaunay/neigh_${name}.dat"
		    params_file="${in_path}/siminfo/params_${name}.txt"
		    lx=$(grep 'lx = ' $params_file | awk '{print $3}')
		    ly=$(grep 'ly = ' $params_file | awk '{print $3}')
		    clx=$(grep 'cellLx = ' $params_file | awk '{print $3}')
		    cly=$(grep 'cellLy = ' $params_file | awk '{print $3}')
		    field_path="${in_path}/cell_field/"
		    field_file="cell-field_${name}_"
		    olap_all_file="${out_path}/local-overlap-all_${name}.dat"
		    olap_neg_file="${out_path}/local-overlap-neg_${name}.dat"
		    olap_pos_file="${out_path}/local-overlap-pos_${name}.dat"
		    olap_log_file="${out_path}/local-overlap-log_${name}.dat"
		    cmd[$jobid]="$olap_exe $N $lx $ly $clx $cly $tstart $tend $tinc $ftinc $pos_file $neigh_file $field_path $field_file $olap_all_file $olap_neg_file $olap_pos_file"
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
