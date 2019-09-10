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
    echo "usage: voronoi.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

vor_py="../src/voronoi.py"

N=100 #100
dt=0.5
Dr=0.0001
obs="hexatic"
tstart=0 # 0
tend=21000000 # 21000000
tinc=10000 # 10000

# Output option
make_movie=1
print_to_screen=0

# Parallel run options
max_jobs=8
cmd=()
jobid=0

# Functions to calculate observables
data_col=0
data_min=0.0
data_max=0.0
data_file=""
remove_file=0 # Remove data file or not after computation

get_asphere() {
    echo "Calculating asphere values ..."
    name=$1
    in_path=$2
    out_path=$3
    data_col=0
    data_min=0.999
    data_max=1.05
    remove_file=1
    shape_file="${in_path}/shape/shape_${name}.dat"
    data_file="${out_path}/asphere-cell_${name}.dat"
    > $data_file # Clear the file
    PI=$(python -c "import math; print '{:.10f}'.format(math.pi)")
    nlines=$(python -c "print $N+2")
    awk -v n=$nlines -v pi=$PI '{
i=(NR-1)%n;
if(i>=2){
as = $1*$1/(4.0*pi*$2);
printf("%.10f\n", as);
} else {print}
}' $shape_file > $data_file
}

get_eccent() {
    echo "Calculating eccent values ..."
    name=$1
    in_path=$2
    out_path=$3
    data_col=0
    data_min=0.0
    data_max=0.5
    remove_file=1
    gyr_file="${in_path}/gyration/gyr_${name}.dat"
    data_file="${out_path}/eccent-cell_${name}.dat"
    > $data_file # Clear the file
    nlines=$(python -c "print $N+2")
    awk -v n=$nlines '{
i=(NR-1)%n;
if(i>=2){
gxx = $1
gyy = $2
gxy = $3
b = (gxx+gyy)/2.0
c = (gxx*gyy-gxy*gxy);
dis = b*b-c
lam1 = b+sqrt(dis)
lam2 = b-sqrt(dis)
ec = (lam1-lam2)/(lam1+lam2)
printf("%.10f\n", ec);
} else {print}
}' $gyr_file > $data_file 
}

get_hexatic() {
    echo "Calculating hexatic values ..."
    name=$1
    in_path=$2
    out_path=$3
    data_col=2
    data_min=-1.0
    data_max=1.0
    remove_file=0
    data_file="${in_path}/hexatic/hex_${name}.dat"
}

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
		    if [[ "$obs" == "asphere" ]]; then
			get_asphere $name $in_path $out_path
		    elif [[ "$obs" == "eccent" ]]; then
			get_eccent $name $in_path $out_path
		    elif [[ "$obs" == "hexatic" ]]; then
			get_hexatic $name $in_path $out_path
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
		    echo "Making voronoi diagrams ... "
		    cmd[$jobid]="python $vor_py $N $lx $ly $Dr $dt $data_col $data_min $data_max $tstart $tend $tinc $make_movie $print_to_screen $pos_file $data_file $out_file"
		    if [[ $remove_file == 1 ]]; then 
			rm $data_file
		    fi
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
