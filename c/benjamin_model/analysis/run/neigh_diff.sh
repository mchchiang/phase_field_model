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

N=100 # 100

while (( $(bc <<< "$d < $d_end") ))
do
    #in_path="${in_dir}/d_${d}/neighbour"
    in_path="${in_dir}/d_${d}/neigh_delaunay"
    if [ -d $in_path ]; then
	out_path="${out_dir}/d_${d}/neigh_delaunay"
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
		    awk -v ncells="$N" 'function abs(v) {return v < 0 ? -v : v} BEGIN {a = 0.0; a2 = 0.0; n = 0; time = 0}{if((NR-1)%(ncells+2)==1){time = $2};if((NR-1)%(ncells+2)>=2){a+=abs(NF-6.0);n+=1};if(n == ncells){print time,a/n;n=0;a=0.0}}' $neigh_file > $neighdiff_file
#		    awk 'function abs(v) {return v < 0 ? -v : v} BEGIN {a = 0.0; a2 = 0.0; n = 0; time = 0}{if((NR-1)%102==1){time = $2};if((NR-1)%102>=2){d=abs(NF-6.0);a+=d; a2+=d*d;n+=1};if(n == 100){a=a/n; a2=a2/n; avg=a; stdev=sqrt((a2-a*a)*n/(n-1.0)); stderr=stdev/sqrt(n); print time,avg,stdev,stderr;n=0;a=0.0;a2=0.0}}' $neigh_file > $neighdiff_file
		fi
	    done
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done
