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
    echo "usage: neigh_num.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=100 # 100

while (( $(bc <<< "$d < $d_end") ))
do
    #in_path="${in_dir}/d_${d}/neighbour"
    #in_path="${in_dir}/d_${d}/neigh_delaunay"
    in_path=$in_dir
    if [ -d $in_path ]; then
	#out_path="${out_dir}/d_${d}/neigh_delaunay"
	out_path=$out_dir
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
		    neighnum_file="${out_path}/neighnum_${name}.dat"
		    neighnumnorm_file="${out_path}/neighnum-norm_${name}.dat"
		    awk -v ncells="$N" '{if((NR-1)%(ncells+2)>=2){print NF} else {print}}' $neigh_file > $neighnum_file
		    awk -v ncells="$N" '{if((NR-1)%(ncells+2)>=2){print NF-6} else {print}}' $neigh_file > $neighnumnorm_file
		fi
	    done
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
