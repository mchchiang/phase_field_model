#!/bin/bash
# init_phase_field_model.sh

deform=$1    # Deformability (d = epsilon/alpha)
peclet=$2    # Peclet number (Pe = v/(R*Dr))
run=$3       # Trial number
run_dir=$4   # Run directory

if [ "$#" != 4 ]; then
    echo "usage: deform peclet run run_dir"
    exit 1
fi

# Format input params
deform=$(python -c "print '{:.3f}'.format($deform)")
peclet=$(python -c "print '{:.3f}'.format($peclet)")

# A function for generating random numbers
max_seed=1000000
function get_rand(){
    # Generate a 4-byte random integer using urandom
    rand=$(od -vAn -N4 -tu4 < /dev/urandom)
    echo $rand
}

#function get_rand(){
#    rand=$(python -c "import random, sys; print random.randint(0,$max_seed)")
#    echo $rand
#}

# Set the model parameters
ncells=100
ncell_x=10
ncell_y=10
confine_radius=8.0
init_radius=7.0
ideal_radius=12.0
cell_Lx=41
cell_Ly=41
phi0=2.0
mu=6000.0
epsilon=0.1
rotate_diff=0.0001
relax_rate=0.1

nsteps=20000000 # 20000000
nequil=10000 # 10000
delta_t=0.5 # 0.5
dump_cm_freq=1000 # 1000
dump_bulk_cm_freq=1000 # 1000
dump_gyr_freq=1000 # 1000
dump_field_freq=100000 # 100000
equildump_cm_freq=1000 # 1000
equildump_gyr_freq=10000 # 10000
equildump_field_freq=10000 # 10000
seed=$(get_rand)

# Set alpha based on deformability
alpha=$(python -c "print '{:f}'.format($epsilon/$deform)")
kappa=$(python -c "print '{:f}'.format($alpha*2.0)")

# Set motility based on Peclet number, rotatioal diff, and ideal radius
motility=$(python -c "print '{:f}'.format($peclet*$rotate_diff*$ideal_radius)")

# Set a hexagonal lattice
tmp_cm_file="cm_$seed.tmp"
tmp_shape_file="shape_$seed.tmp"
size=$(python triangle.py $ncell_x $ncell_y $confine_radius $tmp_cm_file)
Lx=$(echo $size | awk '{print $3}')
Ly=$(echo $size | awk '{print $6}')
python cell_shape.py $cell_Lx $cell_Ly $phi0 $tmp_shape_file circle $init_radius

# Create run directory and set file names
sim_name="cell_N_${ncells}_d_${deform}_Pe_${peclet}_run_${run}"
run_dir="${run_dir}/${sim_name}/"

if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
fi

cm_file="cm_${sim_name}.in"
shape_file="shape_${sim_name}.in"
params_file="params_${sim_name}.txt"
#equildump_cm_file="pos-equil_${sim_name}.dat"
#equildump_gyr_file="gyr-equil_${sim_name}.dat"
#equildump_field_file="field-equil_${sim_name}.dat"
dump_cm_file="pos_${sim_name}.dat"
dump_gyr_file="gyr_${sim_name}.dat"
dump_field_file="field_${sim_name}.dat"
dump_bulk_cm_file="pos-bulk_${sim_name}.dat"

# Copy the template file
params_file=${run_dir}/$params_file
cp params_template.txt $params_file
mv $tmp_cm_file ${run_dir}/$cm_file
mv $tmp_shape_file ${run_dir}/$shape_file

# Replace macros in template with input values
sed -i -- "s/PHI0/${phi0}/g" $params_file
sed -i -- "s/ALPHA/${alpha}/g" $params_file
sed -i -- "s/KAPPA/${kappa}/g" $params_file
sed -i -- "s/MU/${mu}/g" $params_file
sed -i -- "s/EPSILON/${epsilon}/g" $params_file
sed -i -- "s/ROTATE_DIFF/${rotate_diff}/g" $params_file
sed -i -- "s/MOTILITY/${motility}/g" $params_file
sed -i -- "s/RELAX_RATE/${relax_rate}/g" $params_file
sed -i -- "s/IDEAL_RADIUS/${ideal_radius}/g" $params_file

sed -i -- "s/CELL_LX/${cell_Lx}/g" $params_file
sed -i -- "s/CELL_LY/${cell_Ly}/g" $params_file
sed -i -- "s/LX/${Lx}/g" $params_file
sed -i -- "s/LY/${Ly}/g" $params_file

sed -i -- "s/NCELLS/${ncells}/g" $params_file
sed -i -- "s/NSTEPS/${nsteps}/g" $params_file
sed -i -- "s/NEQUIL/${nequil}/g" $params_file
sed -i -- "s/DELTA_T/${delta_t}/g" $params_file
sed -i -- "s/SEED/${seed}/g" $params_file

sed -i -- "s/CM_FILE/${cm_file}/g" $params_file
sed -i -- "s/SHAPE_FILE/${shape_file}/g" $params_file

# Set dumps
function add_dump() {
    params=$1; file=$2;
    if [ "$file" ]; then
	echo "$params $file" >> $params_file 
    fi
}
add_dump "dump_cm $equildump_cm_freq 0 equil" $equildump_cm_file
add_dump "dump_gyr $equildump_gyr_freq 0 equil" $equildump_gyr_file
add_dump "dump_field $equildump_field_freq 0 equil" $equildump_field_file
add_dump "dump_cm $dump_cm_freq 0 main" $dump_cm_file
add_dump "dump_gyr $dump_gyr_freq 0 main" $dump_gyr_file
add_dump "dump_field $dump_field_freq 0 main" $dump_field_file
add_dump "dump_bulk_cm $dump_bulk_cm_freq main" $dump_bulk_cm_file
