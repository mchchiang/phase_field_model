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
#cmin=${10}
#cmax=${11}
in_dir=${10}
out_dir=${11}

if [ "$#" != 11 ]; then
    echo "usage: movie.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc cmin cmax in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

N=100 #100
t_start=0
t_end=100000
t_inc=10000
nframes=$(python -c "print ($t_end-$t_start)/$t_inc+1")
framerate=10

while (( $(bc <<< "$d < $d_end") ))
do
    #in_path="${in_dir}/d_${d}/"
    in_path="${in_dir}/"
    if [ -d $in_path ]; then
	#out_path="${out_dir}/d_${d}/movie/"
	out_path="${out_dir}/snapshot/"
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
		data_file="${in_path}/hexatic_delaunay/hexatic_${name}.dat"
		movie_gp="${out_path}/hexatic_${name}.gp"
		movie="${out_path}/hexatic_${name}"
	    
		echo "
set term png size 640,160 enhanced
set xrange [0:1050]
set yrange [0:1]
set xlabel 'D_rt'
#set ylabel '|{/Symbol Y}_6|'
unset key

#set output '${movie}_t_0.png'
#p '${data_file}' u 1:2 w l lw 5 lc rgb '#4169AF'

do for [i=0:$nframes] {
  set output '${movie}_t_'.i.'.png'
  p '${data_file}' u (\$1*0.00005):(\$1<=i*$t_inc?\$4:1/0) w l lw 1 lc rgb '#4169AF'
}
" > $movie_gp
	    
		echo "Plotting image files for the movie ..."
		gnuplot $movie_gp
		
#		echo "Combining image files into a movie ..."
#		ffmpeg -framerate $framerate -start_number 0 -i ${movie}_t_%d.png -pix_fmt yuv420p ${movie}.mp4
		
#		echo "Cleaning up ..."
#		rm ${movie}_t_*.png
#		rm $movie_gp
	    done
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done


