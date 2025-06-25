#!/bin/bash
#===============================================================================
# Molecular Dynamics Simulation Script for Al Coexistence Studies
#===============================================================================
#
# DESCRIPTION:
#   Script for running molecular dynamics simulations to study aluminum solid-liquid
#   coexistence, with thermalization, partial melting, and equilibrium sampling.
#
# USAGE:
#   ./script.sh
#
# WORKFLOW:
#   Thermalization → Partial melting → Re-thermalization → NVE sampling
#
# REQUIREMENTS:
#   - md.x executable and proper environment for MPI
#
#===============================================================================

# Configuration parameters
TOTAL_ATOMS=7064                    # Total number of atoms in system
TEMP_START=970                      # Starting temperature (K)
TEMP_END=970                        # Ending temperature (K)
TEMP_STEP=1                         # Temperature increment (K)
LIQUID_TEMP=4000                    # Liquid reservoir temperature (K)
ATOMIC_MASS=24.31                   # Atomic mass of Al (amu)
MPI_CORES=4                         # Number of MPI processes

# EAM potential 
    EAM="SIG=9.963404
        SIGREF=3.567764
        MM=3.127240
        EPSILON=0.204844
        NN=6.365543
        -----------"
# Function to write base INCAR parameters
write_incar() {
  cat > INCAR << EOF
        SYSTEM=Al, coexistence
        LSHIFT=.T.; R1REF=5.5; RL=7.321
        NWRITE=0
        tipo=eam
        IBRION=0
        NBLOCK=10
        KBLOCK=100
        POTIM=1.0
        POMASS=$ATOMIC_MASS
        ------------
        $EAM
EOF
}

# Function to add dynamic parameters to INCAR
add_parameters() {
    local steps=$1
    local temperature=$2
    local fixed_atoms=$3
    local use_thermostat=$4
    
    cat >> INCAR << EOF
        NSW=$steps
        TEBEG=$temperature
        ${fixed_atoms:+NFIX=$fixed_atoms}
        ${use_thermostat:+LANDERSON=.T.; NANDERSON=300
}
EOF
}

# Function to run a simulation step
run_simulation_step() {
    local steps=$1
    local temperature=$2
    local fixed_atoms=$3
    local use_thermostat=$4

    write_incar
    add_parameters $steps $temperature "$fixed_atoms" "$use_thermostat"
    
    mpirun -np $MPI_CORES ~/usr/md.x

    # Prepare structure for next step
    if [[ -f CONTCAR ]]; then
        sed -n "1,$((TOTAL_ATOMS + 8))p" CONTCAR > POSCAR
    fi
}

#===============================================================================
# Main simulation loop
#===============================================================================
for temp in $(seq $TEMP_START $TEMP_STEP $TEMP_END); do
    # Step 1: Thermalization
    
    run_simulation_step 5000 $temp "" "thermostat"
    cp INCAR INCAR.1
    
    # Step 2: Partial melting (fix half the atoms)
    fixed_atoms=$((TOTAL_ATOMS / 2))
    run_simulation_step 5000 $LIQUID_TEMP $fixed_atoms "thermostat"
    cp INCAR INCAR.2
    
    # Step 3: Re-thermalization at target temperature
    run_simulation_step 4000 $temp $fixed_atoms "thermostat"
    cp INCAR INCAR.3

    # Step 4: Long NVE evolution for equilibrium sampling
    run_simulation_step 3000 $temp "" ""
    cp INCAR INCAR.4
done
