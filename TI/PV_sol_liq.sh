#!/bin/bash
#===============================================================================
# Volume-Pressure Calculation - Solid and Liquid Phases (FIXED)
#===============================================================================

# Configuration
ELEMENT="Mg"
TEMP=970
NVT_MDLEN=1500
MPI_CORES=6
COEXISTENCE_VOL=21.61

echo "# Volume(Å³) Pressure_Solid(GPa) Pressure_Liquid(GPa)" > VP.dat

# EAM parameters
EAM="SIG=9.963404; SIGREF=3.567764; MM=3.127240; EPSILON=0.204844; NN=6.365543"

# Function to run MD for a specific phase
run_phase() {
    local vol=$1
    local phase=$2
    local poscar_file="POSCAR.$phase"
    
    # Copy and modify POSCAR (FIXED: proper sed syntax)
    cp $poscar_file POSCAR
    sed -i "2s/.*/-$vol/" POSCAR
    
    # Create INCAR and run
    cat > INCAR << EOF
SYSTEM=$ELEMENT ${phase}_VP
LSHIFT=.T.; R1REF=5.5; RL=7.321; NWRITE=0; tipo=eam
IBRION=0; NBLOCK=10; KBLOCK=100; POTIM=1.0; POMASS=24.31
NSW=$NVT_MDLEN; TEBEG=$TEMP; LANDERSON=.T.; NANDERSON=300
$EAM
EOF
    
    mpirun -quiet -np $MPI_CORES ~/usr/md.x > md_${phase}_${vol}.log 2>&1
    
    # Calculate pressure (FIXED: only return the pressure value)
    if [[ -f OUTCAR ]]; then
        pressure=$(awk '/plus kinetic/{avg=($3+$4+$5)/3; sum+=avg; count++} 
                   END{if(count>0) printf "%.2f",(sum/count)/10; else print "N/A"}' OUTCAR)
        
        cp OUTCAR OUTCAR_${phase}_${vol}A3
        echo "$pressure"  # Only return pressure value
    else
        echo "N/A"
    fi
}

# Function to run both phases at one volume
run_volume() {
    local vol=$1
    echo "=== Volume: $vol Å³ ===" >&2  # Send to stderr, not stdout
    
    # Run both phases and capture ONLY pressure values
    p_solid=$(run_phase $vol "sol" 2>/dev/null)  # Suppress stderr
    p_liquid=$(run_phase $vol "liq" 2>/dev/null)  # Suppress stderr
    
    # Clean output to stderr
    echo "Solid: $p_solid GPa, Liquid: $p_liquid GPa" >&2
    
    # Save clean data to VP.dat
    echo "$vol $p_solid $p_liquid" >> VP.dat
}

#===============================================================================
# Main execution
#===============================================================================

# Check required files
for file in POSCAR.sol POSCAR.liq; do
    [[ ! -f $file ]] && { echo "Error: $file not found"; exit 1; }
done

# Get atom count from solid POSCAR
TOTAL_ATOMS=$(sed -n '6p' POSCAR.sol | awk '{print $1}')
EXPECTED_VOL=$(echo "$COEXISTENCE_VOL * 0.986 * $TOTAL_ATOMS" | bc -l)

# Volume range (±10%)
vol_min=$(echo "$EXPECTED_VOL * 0.9" | bc -l)
vol_max=$(echo "$EXPECTED_VOL * 1.1" | bc -l)

echo "Solid-Liquid VP scan: $(printf "%.0f" $vol_min) to $(printf "%.0f" $vol_max) Å³"

# Run calculations
for i in {0..10}; do
    vol=$(echo "$vol_min + $i * ($vol_max - $vol_min) / 10" | bc -l)
    vol=$(printf "%.0f" $vol)
    run_volume $vol
done

echo "Results saved in VP.dat"
cat VP.dat
