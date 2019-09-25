# observables.sh
# Functions to calculate specific observables from raw data

set_asphere_params() {
    name=$1
    in_path=$2
    out_path=$3
    shape_file="${in_path}/shape/shape_${name}.dat"
    data_file="${out_path}/asphere-cell_${name}.dat"
    data_col=0
    data_min=0.999
    data_max=1.05
    tic_start=1.0
    tic_end=1.05
    tic_inc=0.01
    remove_file=1
    args="$shape_file $data_file"
}

set_eccent_params() {
    name=$1
    in_path=$2
    out_path=$3
    gyr_file="${in_path}/gyration/gyr_${name}.dat"
    data_file="${out_path}/eccent-cell_${name}.dat"
    > $data_file # Clear the file
    data_col=0
    data_min=0.0
    data_max=0.5
    tic_start=0.0
    tic_end=0.5
    tic_inc=0.1
    remove_file=1
    args="$gyr_file $data_file"
}

set_hexatic_params() {
    name=$1
    in_path=$2
    out_path=$3
    data_file="${in_path}/hexatic/hexatic-cell_${name}.dat"
    data_col=2
    data_min=-1.0
    data_max=1.0
    tic_start=-1.0
    tic_end=1.0
    tic_inc=0.5
    remove_file=0
    args="$data_file"
}

set_neigh_params() {
    name=$1
    in_path=$2
    out_path=$3
    neigh_file="${in_path}/neighbour/neigh_${name}.dat"
    data_file="${out_path}/neigh-cell_${name}.dat"
    data_col=0
    data_min=4
    data_max=8
    tic_start=4
    tic_end=8
    tic_inc=1
    remove_file=1
    args="$neigh_file $data_file"
}

get_asphere() {
    shape_file=$1
    data_file=$2
    > $data_file # Clear the file
    echo "Calculating asphere values ..."
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
    gyr_file=$1
    data_file=$2
    > $data_file # Clear the file
    echo "Calculating eccent values ..."
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
    # Already precomputed
}

get_neigh() {
    echo "Calculating number of nearest neighbours ..."
    neigh_file=$1
    data_file=$2
    > $data_file # Clear the file
    nlines=$(python -c "print $N+2")
    awk -v n=$nlines '{
i=(NR-1)%n;
if(i>=2){
printf("%d\n", NF);
} else {print}
}' $neigh_file > $data_file 
}
