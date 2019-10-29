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
obs=${10} # asphere, eccent, hexatic, neigh
in_dir=${11}
out_dir=${12}

if [ "$#" != 12 ]; then
    echo "usage: voronoi.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc obs in_dir out_dir"
    exit 1
fi

# Load the functions for calculating the observables
. ./observables.sh --source-only

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

vor_py="../src/voronoi.py"

N=100 #100
dt=0.5
Dr=0.0001
tstart=21000000 # 0
tend=21000000 # 21000000
tinc=10000 # 10000

# Output option
make_movie=0
print_to_screen=1

# Parallel run options
max_jobs=3
cmd=()
jobid=0
job_lot=0

while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/"
    if [ -d $in_path ]; then
	out_path="${out_dir}/d_${d}/voronoi"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		echo "Doing d = $d Pe = $pe run = $run"
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
		pos_file="${in_path}/position/pos_${name}.dat"
		if [ -f $pos_file ]; then
		    # Set observable parameters
		    args=""
		    if [[ "$obs" == "asphere" ]]; then
			set_asphere_params $name $in_path $out_path
			cmd[$jobid]="get_asphere $args"
			jobid=$(bc <<< "$jobid + 1")
		    elif [[ "$obs" == "eccent" ]]; then
			set_eccent_params $name $in_path $out_path
			cmd[$jobid]="get_eccent $args"
			jobid=$(bc <<< "$jobid + 1")
		    elif [[ "$obs" == "hexatic" ]]; then
			set_hexatic_params $name $in_path $out_path
			cmd[$jobid]="get_hexatic $args"
			jobid=$(bc <<< "$jobid + 1")
		    elif [[ "$obs" == "neigh" ]]; then
			set_neigh_params $name $in_path $out_path
			cmd[$jobid]="get_neigh $args"
			jobid=$(bc <<< "$jobid + 1")			
		    else
			continue
		    fi
		    params_file="${in_path}/siminfo/params_${name}.txt"
		    lx=$(grep 'lx = ' $params_file | awk '{print $3}')
		    ly=$(grep 'ly = ' $params_file | awk '{print $3}')
		    if [[ $make_movie == 1 ]]; then
			out_file="${out_path}/vor-${obs}_${name}.mp4"
		    else
			out_file="${out_path}/vor-${obs}_${name}_t_${tstart}.pdf"
		    fi
		    if [[ $remove_file == 1 ]]; then 
			job_lot=3
			cmd[$jobid]="python3 $vor_py $N $lx $ly $Dr $dt $data_col $data_min $data_max $tic_start $tic_end $tic_inc $tstart $tend $tinc $make_movie $print_to_screen $pos_file $data_file $out_file"
			jobid=$(bc <<< "$jobid + 1")
			cmd[$jobid]="rm $data_file"
			jobid=$(bc <<< "$jobid + 1")
		    else
			job_lot=2
			cmd[$jobid]="python3 $vor_py $N $lx $ly $Dr $dt $data_col $data_min $data_max $tic_start $tic_end $tic_inc $tstart $tend $tinc $make_movie $print_to_screen $pos_file $data_file $out_file"
			jobid=$(bc <<< "$jobid + 1")
		    fi
		fi
	    done
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done


# Parallel runs

total_jobs=$jobid

k=0
for (( j=0; j<$job_lot; j++ ))
do
    jobid=$j
    k=0
    while (( $(bc <<< "$jobid < $total_jobs") ))
    do
	for (( i=0; i<$max_jobs && $jobid < $total_jobs; i++ ))
	do
	    jobid=$(bc <<< "$k * $job_lot + $j")
	    echo "${cmd[jobid]} &"
	    ${cmd[jobid]} &
	    k=$(bc <<< "$k + 1")
	    jobid=$(bc <<< "$jobid + 1")
	done
	wait
    done
done
