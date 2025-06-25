#!/bin/bash
#===============================================================================
# Molecular Dynamics Simulation Script for Mg Coexistence Studies
# Advanced Binary Search Algorithm for Efficient Coexistence Detection
#===============================================================================

# Configuration parameters
ELEMENT="Mg"
TOTAL_ATOMS=7064
TEMP_MIN=0                       # Minimum search temperature (K)
TEMP_MAX=970                       # Maximum search temperature (K)
LIQUID_TEMP=4000
ATOMIC_MASS=24.31
NVE_MDLEN=1500                     # Minimum NSW = 1500
MPI_CORES=6

# Search parameters
TEMP_TOLERANCE=10                   # Temperature tolerance for convergence
MAX_ITERATIONS=15                   # Maximum search iterations
REFINEMENT_STEPS=5                  # Refinement steps after finding coexistence

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
    
    # Run MD simulation
    mpirun -quiet -np $MPI_CORES ~/usr/md.x > md_output_${temperature}.log 2>&1

    # Check if simulation completed successfully
    if [[ ! -f CONTCAR ]]; then
        echo "ERROR: CONTCAR not generated for temperature $temperature K"
        return 1
    fi

    # Prepare structure for next step
    sed -n "1,$((TOTAL_ATOMS + 8))p" CONTCAR > POSCAR
    return 0
}

# Function to evaluate a single temperature point
evaluate_temperature() {
    local temp=$1
    local use_long_run=${2:-false}
    
    echo "=== Evaluating temperature: $temp K ===" >&2
    
    # Set simulation lengths (ensure minimum NSW = 1500)
    local therm_steps=2000
    local nve_steps=2000
    if [[ "$use_long_run" == "true" ]]; then
        therm_steps=5000
        nve_steps=3000
        echo "Using long evaluation mode for refinement" >&2
    fi
    
    # FIXED: Always start from fresh original structure
    if [[ -f POSCAR.org ]]; then
        cp POSCAR.org POSCAR
        echo "Reset to original structure" >&2
    else
        echo "ERROR: Original POSCAR not found!" >&2
        return 1
    fi
    
    # Step 1: Thermalization
    echo "Step 1: Thermalization ($therm_steps steps at $temp K)" >&2
    if ! run_simulation_step $therm_steps $temp "" "thermostat"; then
        echo "ERROR: Step 1 failed" >&2
        return 1
    fi
    
    # Step 2: Partial melting (fix half atoms, heat to liquid temp)
    echo "Step 2: Partial melting (fix $((TOTAL_ATOMS / 2)) atoms at $LIQUID_TEMP K)" >&2
    fixed_atoms=$((TOTAL_ATOMS / 2))
    if ! run_simulation_step $therm_steps $LIQUID_TEMP $fixed_atoms "thermostat"; then
        echo "ERROR: Step 2 failed" >&2
        return 1
    fi
    
    # Step 3: Re-thermalization (bring back to test temperature with fixed atoms)
    echo "Step 3: Re-thermalization ($therm_steps steps at $temp K with fixed atoms)" >&2
    if ! run_simulation_step $therm_steps $temp $fixed_atoms "thermostat"; then
        echo "ERROR: Step 3 failed" >&2
        return 1
    fi

    # Step 4: NVE sampling (unfix atoms, equilibrate at test temperature)
    echo "Step 4: NVE sampling ($nve_steps steps at $temp K, no thermostat)" >&2
    if ! run_simulation_step $nve_steps $temp "" ""; then
        echo "ERROR: Step 4 failed" >&2
        return 1
    fi
    
    # Analyze results from Step 4
    density_profile CONTCAR "density_${temp}K.dat"
    local state=$(check_coexistence "density_${temp}K.dat")
    
    echo "Result at $temp K: $state" >&2
    echo "==================================" >&2
    
    # Save results
    cp CONTCAR "CONTCAR_${temp}K"
    cp INCAR "INCAR_${temp}K"
    
    # Return only the state (to stdout)
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
        
        # Evaluate temperature with standard run length  
        local state=$(evaluate_temperature $mid "false")
        
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to evaluate temperature $mid K"
            return 1
        fi
        
        echo "Final result: $state"
        
        # FIXED: Correct binary search logic based on phase
        if [[ "$state" == "Coexistence" ]]; then
            echo "SUCCESS: Found coexistence at $mid K!"
            COEXISTENCE_TEMP=$mid
            SEARCH_SUCCESS=1
            
            # Perform refinement with longer runs
            echo "Performing refinement with longer simulation..."
            local refined_state=$(evaluate_temperature $mid "true")
            echo "Refined result: $refined_state"
            
            if [[ "$refined_state" == "Coexistence" ]]; then
                echo "Coexistence confirmed with longer simulation - SEARCH COMPLETE!"
                return 0
            else
                echo "Refinement showed different phase: $refined_state"
                # Continue search based on refined result
                if [[ "$refined_state" == "Solid-only" ]]; then
                    echo "Refined: Still solid -> need HIGHER temperature"
                    low=$mid
                else  # Liquid-only
                    echo "Refined: Now liquid -> need LOWER temperature"
                    high=$mid
                fi
            fi
        elif [[ "$state" == "Solid-only" ]]; then
            echo "SOLID phase detected -> need HIGHER temperature (melting point is above $mid K)"
            low=$mid  # This moves search range UP: next range will be ($mid to $high)
        else  # Liquid-only
            echo "LIQUID phase detected -> need LOWER temperature (melting point is below $mid K)"  
            high=$mid  # This moves search range DOWN: next range will be ($low to $mid)
        fi
        
        # Check convergence
        local range=$((high - low))
        echo "Current search range: $range K"
        echo "NEW search range for next iteration: $low-$high K"
        
        if [[ $range -le $TEMP_TOLERANCE ]]; then
            echo "Converged to within tolerance ($TEMP_TOLERANCE K)"
            # Final check at estimated coexistence temperature
            local final_temp=$(( (low + high) / 2 ))
            echo "Final verification at $final_temp K with long simulation..."
            local final_state=$(evaluate_temperature $final_temp "true")
            
            if [[ "$final_state" == "Coexistence" ]]; then
                COEXISTENCE_TEMP=$final_temp
                SEARCH_SUCCESS=1
                return 0
            else
                echo "Final verification: $final_state (not coexistence)"
                COEXISTENCE_TEMP=$final_temp
                SEARCH_SUCCESS=0
                return 1
            fi
        fi
        
        iteration=$((iteration + 1))
        echo "-------------------------------------------------"
    done
    
    # If maximum iterations reached
    echo "Maximum iterations reached without finding exact coexistence"
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
echo "Minimum NSW: $NVE_MDLEN"
echo ""

# Backup original structure
if [[ -f POSCAR ]]; then
    cp POSCAR POSCAR.org
    echo "Original POSCAR backed up as POSCAR.org"
else
    echo "ERROR: No POSCAR file found!"
    exit 1
fi

# Perform the search
if binary_search_coexistence; then
    echo ""
    echo "==============================================="
    echo "SEARCH COMPLETED SUCCESSFULLY"
    echo "==============================================="
    echo "Coexistence temperature found: $COEXISTENCE_TEMP K"
    echo "Final structure saved as: CONTCAR_${COEXISTENCE_TEMP}K"
    echo "Final density profile: density_${COEXISTENCE_TEMP}K.dat"
else
    echo ""
    echo "==============================================="
    echo "SEARCH COMPLETED WITH ESTIMATE"
    echo "==============================================="
    echo "Coexistence temperature estimate: $COEXISTENCE_TEMP K"
    echo "Note: Exact coexistence not found within tolerance"
    echo "Consider adjusting search parameters or tolerance"
fi

echo ""
echo "All intermediate results saved with temperature labels"
echo "Check md_output_*K.log files for detailed simulation logs"