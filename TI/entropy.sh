#!/bin/bash
#===============================================================================
# Volume-Pressure Calculation - Manipulates existing POSCAR
#===============================================================================

# Configuration
ELEMENT="Mg"
TEMP=998
NVT_MDLEN=1500
MPI_CORES=6
COEXISTENCE_VOL=23.61
TOTAL_ATOMS=$(sed -n '6p' POSCAR | awk '{print $1}')
EXPECTED_VOL=$(echo "$COEXISTENCE_VOL * 0.986 * $TOTAL_ATOMS" | bc -l)

echo "# Volume(Å³) Pressure(GPa)" > VP.dat

# EAM parameters
EAM="SIG=9.963404; SIGREF=3.567764; MM=3.127240; EPSILON=0.204844; NN=6.365543"

# Function to set volume and run MD
run_volume() {
    local vol=$1
    echo "=== Volume: $vol Å³ ==="
    
    # Modify POSCAR volume (line 2)
    sed -i "2s/.*/-$vol/" POSCAR
    
    # Create INCAR and run
    cat > INCAR << EOF
          SYSTEM=$ELEMENT VP_calc
          LSHIFT=.T.; R1REF=5.5; RL=7.321; NWRITE=0; tipo=eam
          IBRION=0; NBLOCK=10; KBLOCK=100; POTIM=1.0; POMASS=24.31
          NSW=$NVT_MDLEN; TEBEG=$TEMP; LANDERSON=.T.; NANDERSON=300
          $EAM
EOF
    
    mpirun -np $MPI_CORES ~/usr/md.x 
    # Calculate and save pressure
    if [[ -f OUTCAR ]]; then
        pressure=$(awk '/plus kinetic/{avg=($3+$4+$5)/3; sum+=avg; count++} 
                   END{if(count>0) printf "%.2f",(sum/count)/10; else print "N/A"}' OUTCAR)
        
        if [[ "$pressure" != "N/A" ]]; then
            echo "$vol $pressure" >> VP.dat
            echo "Pressure: $pressure GPa"
        fi
        
        cp OUTCAR OUTCAR_${vol}A3
    fi
}

#===============================================================================
# Main execution
#===============================================================================

# Backup original POSCAR
cp POSCAR POSCAR.backup

# Volume range (±10%)
vol_min=$(echo "$EXPECTED_VOL * 0.9" | bc -l)
vol_max=$(echo "$EXPECTED_VOL * 1.1" | bc -l)

echo "Volume-Pressure scan: $(printf "%.0f" $vol_min) to $(printf "%.0f" $vol_max) Å³"

# Run calculations
for i in {0..10}; do
    vol=$(echo "$vol_min + $i * ($vol_max - $vol_min) / 10" | bc -l)
    vol=$(printf "%.0f" $vol)
    run_volume $vol
done

echo "Results saved in VP.dat"
cat VP.dat
