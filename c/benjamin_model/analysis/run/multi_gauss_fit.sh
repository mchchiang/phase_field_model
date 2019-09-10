#!/bin/bash
d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
tstart=$7 #1000000
tend=$8 #21000000
tinc=$9 #1000
run=${10} #1
obs=${11}  #"asphere"
in_dir=${12}
out_dir=${13}

if [ "$#" != 13 ]; then
    echo "usage: multi_gauss_fit.sh d_start d_end d_inc pe_start pe_end pe_inc tstart tend tinc run obs in_dir out_dir"
    exit 1
fi

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

d=$(python -c "print '%.3f' % ($d_start)")
pe=$(python -c "print '%.3f' % ($pe_start)")

guass_R="../src/multi_gauss_fit.R"

N=100 # 100

#overall_gauss_file="${out_dir}/${obs}-modality_cell_N_${N}_d_${d_start}-${d_end}_Pe_${pe_start}-${pe_end}_t_${tstart}-${tend}_run_${run}.dat"
overall_gauss_file="${out_dir}/${obs}-modality_cell_N_${N}_d_${d_start}-${d_end}_Pe_${pe_start}-${pe_end}_t_${tstart}-${tend}.dat"
> $overall_gauss_file

while (( $(bc <<< "$d < $d_end") ))
do
    in_path="${in_dir}/d_${d}/"
    if [ -d $in_path ]; then
	out_path="${out_dir}"
	if [ ! -d $out_path ]; then
	    mkdir -p $out_path
	fi
	pe=$(python -c "print '%.3f' % ($pe_start)")
	while (( $(bc <<< "$pe < $pe_end") ))
	do
#	    name="cell_N_${N}_d_${d}_Pe_${pe}_run_${run}"
	    name="cell_N_${N}_d_${d}_Pe_${pe}"
#	    asphere_file="${in_path}/${obs}/${obs}_${name}.dat"
	    asphere_file=${in_path}/${obs}/${obs}_${name}_run_*.dat
	    asphere_file_1=${in_path}/${obs}/${obs}_${name}_run_1.dat
	    if [ -f $asphere_file_1 ]; then
		echo "d = $d Pe = $pe run = $run"
		gauss_file="${out_path}/${obs}-distrb-gauss-fit_${name}_t_${tstart}-${tend}.dat"
		time_col=1
		val_col=2
		ncomp=2 # Number of modes (number of Gaussians)
		Rscript "../src/multi_gauss_fit.R" $time_col $val_col $tstart $tend $asphere_file > $gauss_file
		means=$(grep "mean:" $gauss_file | awk '{print $2,$3}')
		stdevs=$(grep "stdev:" $gauss_file | awk '{print $2,$3}')
		weights=$(grep "weight:" $gauss_file | awk '{print $2,$3}')
		modality=$(grep "modality:" $gauss_file | awk '{if ($2=="unimodal"){$2=0} else {$2=1}; print $2}')
#		echo "mean: $means"
#		echo "stdev: $stdevs"
#		echo "weight: $weights"
#		echo "modality: $modality"
		echo "$d $pe $means $stdevs $weights $modality" >> $overall_gauss_file
		rm $gauss_file
	    fi
	    pe=$(python -c "print '%.3f' % ($pe + $pe_inc)")
	done
	echo "" >> $overall_gauss_file
    fi
    d=$(python -c "print '%.3f' % ($d + $d_inc)")
done

