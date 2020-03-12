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
    echo "usage: plot_cell_field.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

plot_py="../src/plot_cell_field.py"

N=100 #100
#tstart=12000000
#tend=12000000
tstart=0 #5000000 #20000000
tend=21000000 #5000000 #20000000
tinc=10000
make_movie=1
print_to_screen=0
dr=0.0001
dt=0.5
obs="movie" #"local-overlap"
obs_dir="movie" #"local_overlap"
#obs_name="movie" #"local-overlap"
#obs="local-overlap"
#obs_dir="local_overlap"
obs_name=""
data_col="-1"
data_min=0
data_max=100
tic_start=0
tic_end=100
tic_inc=20

max_jobs=1 # 8
cmd=()
jobid=0

while (( $(bc <<< "$d < $d_end") ))
do
    #in_path="${in_dir}/d_${d}/"
    in_path="${in_dir}"
    if [ -d $in_path ]; then
	#out_path="${out_dir}/d_${d}/snapshot/"
	out_path="${out_dir}/snapshot/"
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
		    params_file="${in_path}/siminfo/params_${name}.txt"
		    lx=$(grep 'lx = ' $params_file | awk '{print $3}')
		    ly=$(grep 'ly = ' $params_file | awk '{print $3}')
		    clx=$(grep 'cellLx = ' $params_file | awk '{print $3}')
		    cly=$(grep 'cellLy = ' $params_file | awk '{print $3}')
		    field_path="${in_path}/cell_field/"
		    field_file="cell-field_${name}_"
		    data_file="${in_path}/${obs_dir}/${obs}_${name}.dat"
		    #out_file="${out_path}/snapshot-field-${obs_name}_${name}_t_${tstart}.pdf"
		    #out_file="${out_path}/snapshot-field_${name}_t_${tstart}.pdf"
		    out_file="${out_path}/movie_${name}.mp4"
		    cmd[$jobid]="python3 $plot_py $N $lx $ly $clx $cly $dr $dt $data_col $data_min $data_max $tic_start $tic_end $tic_inc $tstart $tend $tinc $make_movie $print_to_screen $field_path $field_file $pos_file $data_file $out_file"
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
