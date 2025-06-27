#!/bin/bash

# Sampling Integration Script - Fixed Version
NUM_ATOMS=432
CORES=4
ELEMENT="Mg"

# Create POSCAR.head from existing POSCAR - FIXED
create_poscar_head() {
    if [[ "$1" == "DFT" ]]; then
        # For DFT: Take first 5 lines (comment + volume + lattice), then add element, atoms, Direct
        { head -5 "../$2"; echo "  $ELEMENT"; echo "  $NUM_ATOMS"; echo " Direct"; } > POSCAR.head
    else
        # For EAM: Take first 7 lines (comment + volume + lattice + atoms + Direct)
        head -7 "../$2" > POSCAR.head
    fi
}

# Setup directories and files
setup_sampling() {
    mkdir -p s${NUM_ATOMS}DFT l${NUM_ATOMS}DFT l${NUM_ATOMS}EAM s${NUM_ATOMS}EAM
    [[ -f XDATCAR.sol ]] && sed -i '1,5d' XDATCAR.sol
    [[ -f XDATCAR.liq ]] && sed -i '1,5d' XDATCAR.liq
    for dir in s${NUM_ATOMS}EAM l${NUM_ATOMS}EAM; do cp INCAR $dir/ && sed -i 's/NSW=.*/NSW=0/' $dir/INCAR; done
    for dir in s${NUM_ATOMS}DFT l${NUM_ATOMS}DFT; do cp INCAR $dir/; done
}

# Run integration
run_integration() {
    create_poscar_head "$2" "$3"
    imax=$(wc -l ../$1 | awk -v atoms=$((NUM_ATOMS+1)) '{print int($1/atoms)}')
    
    for ((i=1; i<=imax; i++)); do
        [[ -f OUTCAR.$i ]] && continue
        cp POSCAR.head POSCAR
        awk -v i=$i -v atoms=$((NUM_ATOMS+1)) \
            '{if(NR>(1+(i-1)*atoms) && (NR<=((i-1)*atoms+atoms))) print $1,$2,$3}' \
            ../$1 >> POSCAR
        
        [[ "$2" == "DFT" ]] && mpirun -np $CORES vasp_std > out.log 2>&1 || mpirun -np $CORES md.x > out.log 2>&1
        cp OUTCAR OUTCAR.$i
    done
}

# Main execution
setup_sampling
cd s${NUM_ATOMS}DFT && run_integration "XDATCAR.sol" "DFT" "POSCAR.sol" && cd ..
cd l${NUM_ATOMS}DFT && run_integration "XDATCAR.liq" "DFT" "POSCAR.liq" && cd ..
cd s${NUM_ATOMS}EAM && run_integration "XDATCAR.sol" "EAM" "POSCAR.sol" && cd ..
cd l${NUM_ATOMS}EAM && run_integration "XDATCAR.liq" "EAM" "POSCAR.liq" && cd ..
echo "Done!"

