#!/bin/bash
#===============================================================================
# Volume-Pressure Calculation - Solid and Liquid Phases (with melting)
#===============================================================================

# Configuration
ELEMENT="Mg"
TEMP=970
MELT_TEMP=4000                      # High temperature for melting
MELT_STEPS=4000                     # Steps for melting
NVT_MDLEN=30000
MPI_CORES=6
COEXISTENCE_VOL=21.61

echo "# Volume(Å³) Pressure_Solid(GPa) Pressure_Liquid(GPa)" > VP.dat

# EAM parameters
EAM=$(head -n 1 ../Fit_ls/mdincar)

# Function to melt crystal structure and create liquid POSCAR
create_liquid_poscar() {
    local vol=$1
    echo "Creating liquid structure at volume $vol Å³..." >&2
    
    # Start with crystal structure
    cp POSCAR.sol POSCAR
    sed -i "2s/.*/-$vol/" POSCAR
    
    # Create high-temperature INCAR for melting
    cat > INCAR << EOF
        SYSTEM=$ELEMENT melting
        LSHIFT=.T.; R1REF=5.5; RL=7.321; NWRITE=0; tipo=eam
        IBRION=0; NBLOCK=10; KBLOCK=100; POTIM=1.0; POMASS=24.31
        NSW=$MELT_STEPS; TEBEG=$MELT_TEMP; LANDERSON=.T.; NANDERSON=300
        $EAM
EOF
    
    # Run melting simulation
    mpirun -quiet -np $MPI_CORES ~/usr/md.x > /dev/null 2>&1
    
    # Prepare melted structure for liquid calculation
    if [[ -f CONTCAR ]]; then
        sed -n "1,$((TOTAL_ATOMS + 8))p" CONTCAR > POSCAR
        cp POSCAR POSCAR.liq
        echo "Liquid structure created" >&2
        return 0
    else
        echo "ERROR: Melting failed" >&2
        return 1
    fi
}

# Function to run MD for a specific phase
run_phase() {
    local vol=$1
    local phase=$2
    
    if [[ "$phase" == "sol" ]]; then
        # Use crystal structure
        cp POSCAR.sol POSCAR
        sed -i "2s/.*/-$vol/" POSCAR
    else
        # Create and use liquid structure
        if ! create_liquid_poscar $vol; then
            echo "N/A"
            return 1
        fi
    fi
    
    # Create INCAR for production run
    cat > INCAR << EOF
          SYSTEM=$ELEMENT ${phase}_VP
          LSHIFT=.T.; R1REF=5.5; RL=7.321; NWRITE=0; tipo=eam
          IBRION=0; NBLOCK=10; KBLOCK=100; POTIM=1.0; POMASS=24.31
          NSW=$NVT_MDLEN; TEBEG=$TEMP; LANDERSON=.T.; NANDERSON=300
          $EAM
EOF
    
    # Run production simulation
    mpirun -quiet -np $MPI_CORES ~/usr/md.x > /dev/null 2>&1
    
    # Calculate pressure
    if [[ -f OUTCAR ]]; then
        pressure=$(awk '/plus kinetic/{avg=($3+$4+$5)/3; sum+=avg; count++} 
                   END{if(count>0) printf "%.2f",(sum/count)/10; else print "N/A"}' OUTCAR)
        
        cp OUTCAR 
        echo "$pressure"
    else
        echo "N/A"
    fi
}

# Function to run both phases at one volume
run_volume() {
    local vol=$1
    echo "=== Volume: $vol Å³ ===" >&2
    
    # Run solid phase
    echo "Running solid phase..." >&2
    p_solid=$(run_phase $vol "sol" 2>/dev/null)
    
    # Run liquid phase (includes melting step)
    echo "Running liquid phase (with melting at $MELT_TEMP K)..." >&2
    p_liquid=$(run_phase $vol "liq" 2>/dev/null)
    
    echo "Solid: $p_solid GPa, Liquid: $p_liquid GPa" >&2
    echo "$vol $p_solid $p_liquid" >> VP.dat
}

#===============================================================================
# Main execution
#===============================================================================

# Check required files
cp ../POSCAR POSCAR.sol
[[ ! -f POSCAR.sol ]] && { echo "Error: POSCAR.sol not found"; exit 1; }

# Get atom count from solid POSCAR
TOTAL_ATOMS=$(sed -n '6p' POSCAR.sol | awk '{print $1}')
EXPECTED_VOL=$(echo "$COEXISTENCE_VOL * 0.986 * $TOTAL_ATOMS" | bc -l)

# Volume range (±10%)
vol_min=$(echo "$EXPECTED_VOL * 0.9" | bc -l)
vol_max=$(echo "$EXPECTED_VOL * 1.1" | bc -l)

echo "Solid-Liquid VP scan: $(printf "%.0f" $vol_min) to $(printf "%.0f" $vol_max) Å³"
echo "Melting temperature: $MELT_TEMP K for $MELT_STEPS steps"

# Run calculations
for i in {0..10}; do
    vol=$(echo "$vol_min + $i * ($vol_max - $vol_min) / 10" | bc -l)
    vol=$(printf "%.0f" $vol)
    run_volume $vol
done

echo "Results saved in VP.dat"
cat VP.dat

