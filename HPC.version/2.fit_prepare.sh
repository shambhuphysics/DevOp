#!/bin/bash
#===============================================================================
# VASP AIMD Fitting Input Preparation
#===============================================================================

prepare_fitting_input() {
    echo "Preparing fitting input files..."
    
    # Check if required directories exist
    if [[ ! -d "Solid" ]] || [[ ! -d "Liquid" ]]; then
        echo "Error: Solid and/or Liquid directories not found"
        exit 1
    fi
    
    # Create fitting directory
    mkdir -p Fit_ls
    cd Fit_ls
    
    # Copy and modify XDATCAR files
    cp ../Solid/XDATCAR XDATCAR.s && cp ../Solid/POSCAR .
    sed -i '1,2d' XDATCAR.s
    
    cp ../Liquid/XDATCAR XDATCAR.l
    sed -i '1,7d' XDATCAR.l
    
    # Combine XDATCAR files
    cat XDATCAR.s XDATCAR.l > XDATCAR
    
    # Combine free energy and pressure files
    cat ../Solid/free_energy ../Liquid/free_energy > free_energy
    cat ../Solid/pressure ../Liquid/pressure > pressure
    
    # Count lines in free_energy file
    num_of_conf=$(wc -l free_energy | awk '{print $1}')
    
    # Extract number of atoms from XDATCAR
    num_atoms=$(awk '/^[[:space:]]*[A-Z][a-z]*[[:space:]]*$/ {getline; print $1; exit}' XDATCAR)
    
    # Create ltemp1 file
    cat > ltemp1 << EOF
$num_atoms   0 $num_of_conf   1 15000  5.50  0.0010
0.146646  3.334559  7.115961  4.712721 10.623535
0.000000  0.000000  0.000000  0.000000
0.000000  0.000000  0.000000  0.000000  0.000000
0.50000  0.50000  0.50000  0.50000  5.000000
0.000000  0.000000  0.000000  0.000000
0.000000  0.000000  0.000000  0.000000  0.000000
CHISQ =   0.5518E-03  0.7556E-01  eV^2/at  Iter =  2   14140
EOF
    
    cd ..
    echo "Fitting input files prepared in Fit_ls/ directory"
    echo "Number of atoms: $num_atoms, Number of configurations: $num_of_conf"
}

# Run the function
prepare_fitting_input
