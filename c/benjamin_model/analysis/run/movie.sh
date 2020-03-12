#!/bin/bash
# movie.sh
# Combine field output in different time steps into a movie

d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
run_start=$7
run_end=$8
run_inc=$9
cmin=${10}
cmax=${11}
in_dir=${12}
out_dir=${13}

if [ "$#" != 13 ]; then
    echo "usage: movie.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc cmin cmax in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=36 #100
t_start=0
t_end=1000000 # 21000000
t_inc=1000 # 100000
zeros=${t_inc:1}
nframes=$(python -c "print ($t_end-$t_start)/$t_inc+1")
framerate=5

while (( $(bc <<< "$d < $d_end") ))
do
#    in_path="${in_dir}/d_${d}/"
    in_path="${in_dir}/"
    if [ -d $in_path ]; then
#	out_path="${out_dir}/d_${d}/movie/"
	out_path="${out_dir}/"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		echo "Doing N = ${N} d = ${d} Pe = ${pe}"
		name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
#		field_file="${in_path}/field/field_${name}.dat"
#		field_file="${in_path}/field_tmp/field_${name}.dat"
		field_file="${in_path}/field_${name}.dat"
#		param_file="${in_path}/siminfo/params_${name}.txt"
		param_file="${in_path}/params_${name}.txt"
	    
		lx=$(grep 'lx = ' $param_file | awk '{print $3}')
		ly=$(grep 'ly = ' $param_file | awk '{print $3}')
		lxm1=$(python -c "print $lx-1")
		lym1=$(python -c "print $ly-1")
	    
		movie_gp="${out_path}/movie_${name}.gp"
		movie="${out_path}/movie_${name}"
	    
		echo "
set term png size 800,600
set xrange [0:$lxm1]
set yrange [0:$lym1]
set cbrange [$cmin:$cmax]
unset key

set size ratio -1
set output '${movie}_t_0.png'
set label 1 't = 0' front at graph 0.8, 0.95
p '${field_file}.0' u 1:2:3 w image

do for [i=1:$nframes] {
  set size ratio -1
  set output '${movie}_t_'.i.'.png'
  set label 1 't = '.i front at graph 0.8, 0.95
  p '${field_file}.'.i.'$zeros' u 1:2:3 w image
}
" > $movie_gp
	    
		echo "Plotting image files for the movie ..."
		gnuplot $movie_gp
		
		echo "Combining image files into a movie ..."
		ffmpeg -framerate $framerate -start_number 0 -i ${movie}_t_%d.png -pix_fmt yuv420p ${movie}.mp4
		
		echo "Cleaning up ..."
		rm ${movie}_t_*.png
		rm $movie_gp
	    done
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done


