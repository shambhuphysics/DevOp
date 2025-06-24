#!/bin/bash
#===============================================================================
# VASP AIMD Melting Workflow
#===============================================================================

# Configuration
ELEMENT="${ELEMENT:-Al}"
ATOMS="${ATOMS:-4}"
KPOINTS_MESH="${KPOINTS_MESH:-2 2 2}"
MDSTEPS_L="${MDSTEPS_L:-500}"
MDSTEPS_M="${MDSTEPS_M:-1000}"
TM="${TM:-920}"
TL="${TL:-4000}"
NPROC="${NPROC:-4}"

# Calculate electronic temperatures
TM_ELEC=$(echo "scale=6; $TM * 8.617e-5" | bc)
TL_ELEC=$(echo "scale=6; $TL * 8.617e-5" | bc)

# Functions
write_incar() {
    cat > INCAR << EOF
SYSTEM = $ELEMENT, $ATOMS atoms
MAXMIX = 60
NPAR = 4
LCHARG = .FALSE.
LWAVE = .FALSE.
NWRITE = 0
IALGO = 48
ISYM = 0
IBRION = 0
NBLOCK = 5
KBLOCK = 1
SMASS = 1.0
POTIM = 2.0
ISIF = 2
ISMEAR = -1
EOF
}

add_temp() {
    echo "TEBEG = $1" >> INCAR
    echo "SIGMA = $2" >> INCAR
    echo "NSW = $3" >> INCAR
}

create_kpoints() {
    cat > KPOINTS << EOF
K-Points
0
Monkhorst Pack
$KPOINTS_MESH
0 0 0
EOF
}

copy_potcar() {
    cp ~/POTs/$ELEMENT POTCAR
}

# Main workflow
echo "Starting VASP AIMD melting simulation..."

# Setup liquid phase
mkdir -p Liquid
cp POSCAR Liquid/
cd Liquid

# Generate input files
write_incar
add_temp $TL $TL_ELEC $MDSTEPS_L
create_kpoints
copy_potcar

# Run high-T equilibration
echo "Running liquid equilibration at ${TL}K..."
mpirun -np $NPROC vasp_std

# Prepare melting simulation
cp CONTCAR POSCAR
write_incar
add_temp $TM $TM_ELEC $MDSTEPS_M

# Run melting simulation
echo "Running melting simulation at ${TM}K..."
mpirun -np $NPROC vasp_std

echo "AIMD melting simulation completed!"
