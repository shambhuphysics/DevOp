#!/bin/bash
#===============================================================================
# Molecular Dynamics Simulation Script for Al Coexistence Studies
# Advanced Binary Search Algorithm for Efficient Coexistence Detection
#===============================================================================

# Configuration parameters
ELEMENT="Mg"
TOTAL_ATOMS=7064
TEMP_MIN=1500                        # Minimum search temperature (K)
TEMP_MAX=2000                       # Maximum search temperature (K)
LIQUID_TEMP=4000
ATOMIC_MASS=24.31
NVE_MDLEN=500
MPI_CORES=4

# Search parameters
TEMP_TOLERANCE=10                   # Reduced tolerance for better precision
MAX_ITERATIONS=10                   # Reduced iterations to avoid excessive runs
REFINEMENT_STEPS=5                  # Reduced refinement steps

# Global variables for results
COEXISTENCE_TEMP=0
SEARCH_SUCCESS=0

# EAM potential 
EAM="SIG=9.963404
     SIGREF=3.567764
     MM=3.127240
     EPSILON=0.204844
     NN=6.365543"

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

# Function to monitor solid, liquid or coexistence
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
    
    # Fixed: Complete mpirun command and redirect output to avoid unit conflicts
    mpirun -quiet -np $MPI_CORES ~/usr/md.x > md_output_${temperature}.log 2>&1

    # Prepare structure for next step
    if [[ -f CONTCAR ]]; then
        sed -n "1,$((TOTAL_ATOMS + 8))p" CONTCAR > POSCAR
    fi
}

# Function to evaluate a single temperature point
evaluate_temperature() {
    local temp=$1
    local use_short=${2:-false}
    
    echo "=== Evaluating temperature: $temp K ==="
    
    # Use shorter runs for initial bracketing
    local therm_steps=3000   # Reduced steps to avoid long runs
    local nve_steps=3000
    if [[ "$use_short" == "true" ]]; then
        therm_steps=3000
        nve_steps=2000
        echo "Using short evaluation mode"
    fi
    
    # Reset to original structure for each temperature
    if [[ -f POSCAR.org ]]; then
        cp POSCAR.org POSCAR
    fi
    
    # Step 1: Thermalization
    echo "Step 1: Thermalization ($therm_steps steps)"
    run_simulation_step $therm_steps $temp "" "thermostat"
    
    # Step 2: Partial melting
    echo "Step 2: Partial melting"
    fixed_atoms=$((TOTAL_ATOMS / 2))
    run_simulation_step $therm_steps $LIQUID_TEMP $fixed_atoms "thermostat"
    
    # Step 3: Re-thermalization
    echo "Step 3: Re-thermalization"
    run_simulation_step $therm_steps $temp $fixed_atoms "thermostat"

    # Step 4: NVE sampling
    echo "Step 4: NVE sampling ($nve_steps steps)"
    run_simulation_step $nve_steps $temp "" ""
    
    # Analyze results
    density_profile CONTCAR "density_${temp}K.dat"
    local state=$(check_coexistence "density_${temp}K.dat")
    
    echo "Result at $temp K: $state"
    echo "=================================="
    
    # Save results
    cp CONTCAR "CONTCAR_${temp}K"
    cp INCAR "INCAR_${temp}K"
    
    echo "$state"
}

# Function to perform binary search for coexistence
binary_search_coexistence() {
    local low=$TEMP_MIN
    local high=$TEMP_MAX
    local iteration=1
    
    echo "Starting binary search between $low K and $high K"
    echo "Target tolerance: $TEMP_TOLERANCE K"
    echo "================================================="
    
    while [[ $iteration -le $MAX_ITERATIONS ]]; do
        local mid=$(( (low + high) / 2 ))
        
        echo "Iteration $iteration: Testing mid-point $mid K (range: $low-$high K)"
        
        local state=$(evaluate_temperature $mid "true")
        
        echo "Temperature $mid K -> $state"
        
        # FIXED: Correct binary search logic
        if [[ "$state" == "Coexistence" ]]; then
            echo "SUCCESS: Found coexistence at $mid K!"
            COEXISTENCE_TEMP=$mid
            SEARCH_SUCCESS=1
            return 0
        elif [[ "$state" == "Solid-only" ]]; then
            echo "Solid phase -> searching HIGHER temperatures"
            low=$mid
        else  # Liquid-only
            echo "Liquid phase -> searching LOWER temperatures"
            high=$mid
        fi
        
        # Check convergence
        local range=$((high - low))
        echo "Current search range: $range K"
        
        if [[ $range -le $TEMP_TOLERANCE ]]; then
            echo "Converged to within tolerance ($TEMP_TOLERANCE K)"
            break
        fi
        
        iteration=$((iteration + 1))
        echo "Next iteration will test range: $low-$high K"
        echo "-------------------------------------------------"
    done
    
    # If no exact coexistence found, return best estimate
    COEXISTENCE_TEMP=$(( (low + high) / 2 ))
    SEARCH_SUCCESS=0
    return 1
}

#===============================================================================
# Main execution
#===============================================================================

echo "Advanced Coexistence Temperature Search"
echo "======================================="
echo "Element: $ELEMENT"
echo "Total atoms: $TOTAL_ATOMS"
echo "Search range: $TEMP_MIN - $TEMP_MAX K"
echo "Tolerance: $TEMP_TOLERANCE K"
echo ""

# Backup original structure
if [[ -f POSCAR ]]; then
    cp POSCAR POSCAR.org
    echo "Original POSCAR backed up"
else
    echo "ERROR: No POSCAR file found!"
    exit 1
fi

# Perform the search
binary_search_coexistence

echo ""
echo "==============================================="
echo "FINAL RESULTS"
echo "==============================================="

if [[ $SEARCH_SUCCESS -eq 1 ]]; then
    echo "Coexistence temperature found: $COEXISTENCE_TEMP K"
    echo "Final structure saved as: CONTCAR_${COEXISTENCE_TEMP}K"
    echo "Final density profile: density_${COEXISTENCE_TEMP}K.dat"
else
    echo "Coexistence temperature estimate: $COEXISTENCE_TEMP K"
    echo "Note: Exact coexistence may require further refinement"
fi

echo ""
echo "Search completed successfully!"


# I am surprised, if at 950K, it just found coexistence why it lower temperature, 
# why it just simply does not start with Long refinced coexistence , 

# when search T= 2000-15000
# also are you copoing the fresh POSCAR.coex into POSCAR , after new calculaiton after step4?