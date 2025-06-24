#!/bin/bash
#===============================================================================
# VASP AIMD Melting Workflow - Liquid and Solid Phases
#===============================================================================
#
# DESCRIPTION:
#   Automated workflow for VASP Ab Initio Molecular Dynamics (AIMD) simulations
#   to study melting behavior. Runs two-stage MD calculations for liquid and/or 
#   solid phases at specified temperatures.
#
# USAGE:
#   ./script.sh liquid             # Run liquid phase only
#   ./script.sh solid              # Run solid phase only
#   ./script.sh both               # Run both liquid and solid phases
#
# WORKFLOW:
#   Liquid: High-T equilibration (TL) → Melting temperature (TM)
#   Solid:  Melting temperature (TM) → Melting temperature (TM)
#
# ENVIRONMENT VARIABLES:
#   ELEMENT      - Element symbol (default: Al)
#   ATOMS        - Number of atoms (default: 4)
#   TM           - Melting temperature in K (default: 920)
#   TL           - Liquid temperature in K (default: 4000)
#   MDSTEPS_L    - MD steps for equilibration (default: 500)
#   MDSTEPS_M    - MD steps for melting (default: 1000)
#   NBLOCKs      - NBLOCK value for VASP (default: 5)
#   NPROC        - Number of processors (default: 4)
#
# REQUIREMENTS:
#   - POSCAR file in current directory
#   - POTCAR files in ~/POTs/ELEMENT/
#   - VASP executable (vasp_std) in PATH
#
#===============================================================================

# Configuration
ELEMENT="${ELEMENT:-Al}"
ATOMS="${ATOMS:-4}"
KPOINTS_MESH="${KPOINTS_MESH:-2 2 2}"
MDSTEPS_L="${MDSTEPS_L:-50}"
MDSTEPS_M="${MDSTEPS_M:-100}"
NBLOCKs="${NBLOCKs:-5}"
TM="${TM:-920}"
TL="${TL:-4000}"
NPROC="${NPROC:-4}"

# Get phase from command line argument
PHASE="${1:-both}"

# Calculate electronic temperatures
TM_ELEC=$(awk "BEGIN {printf \"%.6f\", $TM * 8.617e-5}")
TL_ELEC=$(awk "BEGIN {printf \"%.6f\", $TL * 8.617e-5}")

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
      NBLOCK = $NBLOCKs; KBLOCK = 1
      SMASS = 1.0
      POTIM = 2.0
      ISIF = 1
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
    cp ~/POTs/$ELEMENT/POTCAR .
}

extract_data() {
    local nblock=${1:-$NBLOCKs}
    echo "Extracting data every $nblock steps..."
    
    # Extract free energy
    awk -v nblock="$nblock" '/free  en/ {count++; if (count % nblock == 0) print count/nblock, $5}' OUTCAR > free_energy
    
    # Extract pressure
    awk -v nblock="$nblock" '/external pressure/ {count++; if (count % nblock == 0) print count/nblock, $4}' OUTCAR > pressure
    
    echo "Data extracted to free_energy and pressure files"
}

run_phase() {
    local phase_name=$1
    local temp1=$2
    local temp1_elec=$3
    local temp2=$4
    local temp2_elec=$5
    
    echo "Starting VASP AIMD ${phase_name} simulation..."
    
    # Setup phase directory
    mkdir -p $phase_name
    cp POSCAR $phase_name/
    cd $phase_name
    
    # Generate input files
    write_incar
    add_temp $temp1 $temp1_elec $MDSTEPS_L
    cp INCAR INCAR.1
    create_kpoints
    copy_potcar
    
    # Run first stage
    echo "Running ${phase_name} equilibration at ${temp1}K..."
    mpirun -np $NPROC vasp_std
    extract_data
    
    # Prepare second stage
    cp CONTCAR POSCAR
    write_incar
    add_temp $temp2 $temp2_elec $MDSTEPS_M
    
    # Run second stage
    echo "Running ${phase_name} simulation at ${temp2}K..."
    mpirun -np $NPROC vasp_std
    extract_data
    
    cd ..
}

# Main workflow
case $PHASE in
    "liquid")
        run_phase "Liquid" $TL $TL_ELEC $TM $TM_ELEC
        ;;
    "solid")
        run_phase "Solid" $TM $TM_ELEC $TM $TM_ELEC
        ;;
    "both")
        run_phase "Liquid" $TL $TL_ELEC $TM $TM_ELEC
        run_phase "Solid" $TM $TM_ELEC $TM $TM_ELEC
        ;;
    *)
        echo "Usage: $0 {liquid|solid|both}"
        echo "  liquid - Run liquid phase only"
        echo "  solid  - Run solid phase only"
        echo "  both   - Run both phases (default)"
        exit 1
        ;;
esac


