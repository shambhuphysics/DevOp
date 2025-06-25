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
ELEMENT="Mg"
TOTAL_ATOMS=7064                    # Total number of atoms in system
TEMP_START=500                      # Starting temperature (K)
TEMP_END=1200                       # Maximum temperature (K) - add this
TEMP_STEP=30                        # Temperature increment (K)
LIQUID_TEMP=4000                    # Liquid reservoir temperature (K)
ATOMIC_MASS=24.31                   # Atomic mass of Al (amu)
NVE_MDLEN=500
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
        SYSTEM=$ELEMENT, coexistence
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
        ${use_thermostat:+LANDERSON=.T.; NANDERSON=300}
EOF
}

# Function to check planar density
density_profile() {
    local ndim=400
    local infile="${1:-CONTCAR}"
    local outfile="${2:-density.dat}"
    
    awk -v ndim="$ndim" '
    BEGIN {
        dx = 1.0/ndim
        for(j=1; j<=ndim; j++) d[j] = 0
    }
    NR==2 { volume = $1 }
    NR==6 { natoms = $1 }
    NR>7 && NR<=7+natoms {
        c = $3
        for(j=1; j<=ndim; j++) {
            x = j/ndim - dx
            if(c >= x && c < x+dx) d[j]++
        }
    }
    END {
        for(j=1; j<=ndim; j++) {
            print j, d[j]
        }
    }' "$infile" > "$outfile"
}

#Function to monitor solid, liquid or coexistence;
check_coexistence() {
    awk 'BEGIN{z=0;c=0} $2<0.001{z++;c=0} $2>=0.001{c++} END{
        if(z>5 && c>50) print "Coexistence"
        else if(z>5) print "Solid-only" 
        else print "Liquid-only"
    }' "${1:-density.dat}"
}

# Function to run a simulation step
run_simulation_step() {
    local steps=$1
    local temperature=$2
    local fixed_atoms=$3
    local use_thermostat=$4

    write_incar
    add_parameters $steps $temperature "$fixed_atoms" "$use_thermostat"
    mpirun -np $MPI_CORES ~/usr/md.x >out

    # Prepare structure for next step
    if [[ -f CONTCAR ]]; then
        sed -n "1,$((TOTAL_ATOMS + 8))p" CONTCAR > POSCAR
    fi
}

#===============================================================================
# Main simulation loop with adaptive temperature search
#===============================================================================
current_temp=$TEMP_START

while [[ $current_temp -le $TEMP_END ]]; do
    echo "Testing temperature: $current_temp K"
    
    # Step 1: Thermalization
    echo "Step 1: Thermalization"
    run_simulation_step 2000 $current_temp "" "thermostat"
    cp INCAR INCAR.1
    
    # Step 2: Partial melting (fix half the atoms)
    echo "Step 2: Partial melting"
    fixed_atoms=$((TOTAL_ATOMS / 2))
    run_simulation_step 2000 $LIQUID_TEMP $fixed_atoms "thermostat"
    cp INCAR INCAR.2
    
    # Step 3: Re-thermalization at target temperature
    echo "Step 3: Re-thermalization"
    run_simulation_step 2000 $current_temp $fixed_atoms "thermostat"
    cp INCAR INCAR.3

    # Step 4: Long NVE evolution for equilibrium sampling
    echo "Step 4: NVE sampling"
    density_profile CONTCAR
    run_simulation_step $NVE_MDLEN $current_temp "" ""
    density_profile CONTCAR
    
    # Check the coexistence state
    state=$(check_coexistence)
    echo "Current state at $current_temp K: $state"
    
    cp density.dat d_retherm.dat  
    cp INCAR INCAR.4
    
    # Decision logic based on state
    if [[ "$state" == "Coexistence" ]]; then
        echo "SUCCESS: Found coexistence at $current_temp K!"
        break
    elif [[ "$state" == "Solid-only" ]]; then
        echo "Solid phase detected, increasing temperature"
        current_temp=$((current_temp + TEMP_STEP))
        cp POSCAR.org POSCAR
    else  # Liquid-only
        echo "Liquid phase detected, decreasing temperature"
        current_temp=$((current_temp - TEMP_STEP))
        cp POSCAR.org POSCAR
        # Prevent going below reasonable temperature
        if [[ $current_temp -lt 500 ]]; then
            echo "Temperature too low, stopping search"
            break
        fi
    fi
    
    echo "Next temperature will be: $current_temp K"
    echo "----------------------------------------"
done

if [[ "$state" == "Coexistence" ]]; then
    echo "Coexistence found at temperature: $current_temp K"
else
    echo "Coexistence search completed without finding stable coexistence"
fi