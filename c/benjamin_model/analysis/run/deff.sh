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

avg_py="../src/AverageMultiFiles.py"

if [ "$#" != 11 ]; then
    echo "usage: deff.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
d_old=$d
pe=$(python -c "print '%.3f' % ($pe_start)")

N=100 # 100
tstart=5000000
tend=15000000

# Convert to simulation time (each timestep is 0.5 unit of time)
tmin=$(python -c "print $tstart/2.0")
tmax=$(python -c "print $tend/2.0")

fitgp="msd_fit.gp"
fitlog="msd_fit.log"

deff_overall_file="${out_dir}/deff_cell_N_${N}_d_${d_start}-${d_end}-${d_inc}_Pe_${pe_start}-${pe_end}-${pe_inc}_t_${tstart}-${tend}.dat"
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
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
	    name="cell_N_${N}_d_${d}_Pe_${pe}"
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		msd_file="${in_path}/msd/t_1000000-21000000/msd-tavg_${name}_run_${run}.dat"
		
		if [ -f $msd_file ]; then
		    echo "Doing d = $d Pe = $pe run = $run"
		    > $fitlog
		    > $fitgp
		    deff_file="${out_path}/deff_${name}_run_${run}.dat"
		    params_file="${in_path}/siminfo/params_${name}_run_${run}.txt"
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
		    # Compute the effective diffusion coefficient, normalised 
		    # by the diffusion of a self-propelled particle
		    deff=$(python -c "print ${slope}/4.0/(${v}**2/(2.0*${Dr})) if (${v} > 0) else 0.0")

		    echo "$deff" > $deff_file
		fi
	    done

	    # Average the results from different runs
	    if [ -f "${out_path}/deff_${name}_run_1.dat" ]; then
		deff_avg_file="${out_path}/deff_${name}_avg.dat"
		python $avg_py -1 0 -1 -1 $deff_avg_file "${out_path}/deff_${name}"_run_*.dat
		data=$(cat $deff_avg_file)
		echo "$d $pe $data" >> $deff_overall_file
		rm "${out_path}/deff_${name}"_run_*.dat
		rm $deff_avg_file
	    fi
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done

rm $fitgp
rm $fitlog
