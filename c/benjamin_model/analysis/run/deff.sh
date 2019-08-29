#!/bin/bash
d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
in_dir=$7
out_dir=$8

run=2

if [ "$#" != 8 ]; then
    echo "usage: deff.sh d_start d_end d_inc pe_start pe_end pe_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
d_old=$d
pe=$(python -c "print '%.3f' % ($pe_start)")

N=100 # 100
#tshiftend=20000000 # 20000000
tstart=5000000
tend=15000000

# Convert to simulation time (each timestep is 0.5 unit of time)
tmin=$(python -c "print $tstart/2.0")
tmax=$(python -c "print $tend/2.0")

fitgp="msd_fit.gp"
fitlog="msd_fit.log"

deff_overall_file="${out_dir}/deff_cell_N_${N}_d_${d_start}-${d_end}_Pe_${pe_start}_${pe_end}_t_${tstart}-${tend}_run_${run}.dat"
> $deff_overall_file

while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}"    
    if [ -d $in_path ]; then
	if [ $d != $d_old ]; then
	    echo "" >> $deff_overall_file
	    d_old=$d
	fi
	out_path="${out_dir}/d_${d}/deff"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	deff_file="${out_path}/deff_cell_N_${N}_d_${d}_run_${run}.dat"
	> $deff_file
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}" # avg file for more than 1 run
	    msd_file="${in_path}/msd/msd-tavg_${name}.dat"
#	    msd_file="${in_path}/msd/msd-tavg_${name}_ts_${tshiftend}.dat"
	    
	    > $fitlog
	    > $fitgp
	    
	    params_file="${in_path}/siminfo/params_${name}.txt"
	    v=$(grep 'v = ' $params_file | awk '{print $3}')
	    Dr=$(grep 'Dr = ' $params_file | awk '{print $3}')
	    echo "
f(x)=a*x+b
a=1.0
b=1.0
set fit logfile '${fitlog}'
set fit quiet
fit [${tmin}:${tmax}] f(x) '${msd_file}' u (\$1*0.5):(\$2) via a,b" > $fitgp

	    gnuplot $fitgp
	    # Grab the fitted slope and intercept
	    slope=$(grep -E "a +[=] [0-9.e-]* +[+][/][-]" $fitlog | awk '{print $3}')
	    intercept=$(grep -E "b +[=] [0-9.e-]* +[+][/][-]" $fitlog | awk '{print $3}')
	    echo $slope $intercept
	    # Compute the effective diffusion coefficient, normalised by the
	    # diffusion of a self-propelled particle
	    deff=$(python -c "print ${slope}/4.0/(${v}**2/(2.0*${Dr}))")
	    echo "$d $pe $deff" >> $deff_file
	    echo "$d $pe $deff" >> $deff_overall_file
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done

rm $fitgp
rm $fitlog
