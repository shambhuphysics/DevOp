#!/bin/bash

# Birch-Murnaghan EOS Fitting and Bulk Modulus Calculation
# Author: Shell Script Implementation
# Date: $(date)

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -f FILE    Input VP.dat file (default: VP.dat)"
    echo "  -p PRESS   Target pressure for volume calculation (default: 0.5)"
    echo "  -h         Show this help message"
    exit 1
}

# Default values
INPUT_FILE="VP.dat"
TARGET_PRESSURE=0.5

# Parse command line arguments
while getopts "f:p:h" opt; do
    case $opt in
        f) INPUT_FILE="$OPTARG" ;;
        p) TARGET_PRESSURE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file $INPUT_FILE not found!"
    echo "Creating sample VP.dat file..."
    cat > VP.dat << 'EOF'
# Volume(Å³) Pressure_Solid(GPa) Pressure_Liquid(GPa)
8284 15.74 15.41
8468 12.87 12.07
8653 10.72 10.22
8837 8.90 8.13
9021 5.57 5.19
9205 5.61 4.90
9389 3.62 4.08
9573 2.45 2.74
9757 1.77 0.55
9941 0.77 0.63
10125 -0.79 -0.66
EOF
    echo "Sample VP.dat created. Please modify with your data if needed."
fi

echo "=========================================="
echo "Birch-Murnaghan EOS Fitting Script"
echo "=========================================="
echo "Input file: $INPUT_FILE"
echo "Target pressure: $TARGET_PRESSURE GPa"
echo ""

# Create gnuplot script for fitting Birch-Murnaghan equation
create_gnuplot_fit_script() {
    local phase=$1
    local column=$2
    local output_prefix=$3
    
    cat > fit_${phase}.gp << EOF
# Gnuplot script for fitting Birch-Murnaghan EOS for $phase phase
set terminal png enhanced size 800,600
set output '${output_prefix}_fit.png'

# Define the Birch-Murnaghan equation
# P = (3/2) * B0 * ((V0/V)**(7/3) - (V0/V)**(5/3)) * (1 + (3/4) * (B0_prime - 4) * ((V0/V)**(2/3) - 1))
birch(V, V0, B0, B0_prime) = (3.0/2.0) * B0 * ((V0/V)**(7.0/3.0) - (V0/V)**(5.0/3.0)) * (1 + (3.0/4.0) * (B0_prime - 4) * ((V0/V)**(2.0/3.0) - 1))

# Initial parameter guesses
V0 = 11000
B0 = 100
B0_prime = 4.5

# Fit the data
fit birch(x, V0, B0, B0_prime) '$INPUT_FILE' using 1:$column via V0, B0, B0_prime

# Save fitted parameters to file
set print '${output_prefix}_params.dat'
print "# Fitted Parameters for $phase phase"
print "V0 = ", V0
print "B0 = ", B0
print "B0_prime = ", B0_prime
set print

# Plot the data and fit
set xlabel 'Volume (Å³)'
set ylabel 'Pressure (GPa)'
set title 'Birch-Murnaghan EOS Fit - $phase Phase'
set grid

plot '$INPUT_FILE' using 1:$column with points pt 7 ps 1.5 title '$phase Data', \\
     birch(x, V0, B0, B0_prime) with lines lw 2 title 'Fitted Curve'

# Generate fitted V-P data for interpolation
set table '${output_prefix}_fitted.dat'
set samples 1000
plot [8000:10500] birch(x, V0, B0, B0_prime)
unset table

# Reset terminal
set terminal x11
EOF
}

# Function to extract parameters from gnuplot output
extract_parameters() {
    local param_file=$1
    local prefix=$2
    
    V0=$(grep "V0 =" $param_file | awk '{print $3}')
    B0=$(grep "B0 =" $param_file | awk '{print $3}')
    B0_prime=$(grep "B0_prime =" $param_file | awk '{print $3}')
    
    echo "${prefix}_V0=$V0"
    echo "${prefix}_B0=$B0"
    echo "${prefix}_B0_prime=$B0_prime"
}

# Function to find volume for target pressure using bisection method
find_volume_for_pressure() {
    local target_p=$1
    local V0=$2
    local B0=$3
    local B0_prime=$4
    local phase=$5
    
    # Create awk script for bisection method
    cat > find_volume.awk << 'EOF'
function birch_pressure(V, V0, B0, B0_prime) {
    return (3.0/2.0) * B0 * ((V0/V)^(7.0/3.0) - (V0/V)^(5.0/3.0)) * (1 + (3.0/4.0) * (B0_prime - 4) * ((V0/V)^(2.0/3.0) - 1))
}

function abs(x) { return x < 0 ? -x : x }

BEGIN {
    target = TARGET_P
    V0 = V0_VAL
    B0 = B0_VAL
    B0_prime = B0_PRIME_VAL
    
    # Initial bounds
    V_low = 8000
    V_high = 12000
    tolerance = 1e-6
    max_iter = 1000
    
    for (i = 1; i <= max_iter; i++) {
        V_mid = (V_low + V_high) / 2.0
        P_mid = birch_pressure(V_mid, V0, B0, B0_prime)
        
        if (abs(P_mid - target) < tolerance) {
            printf "%.6f\n", V_mid
            exit
        }
        
        P_low = birch_pressure(V_low, V0, B0, B0_prime)
        
        if ((P_low - target) * (P_mid - target) < 0) {
            V_high = V_mid
        } else {
            V_low = V_mid
        }
    }
    
    printf "%.6f\n", V_mid
}
EOF

    awk -v TARGET_P=$target_p -v V0_VAL=$V0 -v B0_VAL=$B0 -v B0_PRIME_VAL=$B0_prime -f find_volume.awk
}

# Function to calculate bulk modulus
calculate_bulk_modulus() {
    local V=$1
    local V0=$2
    local B0=$3
    local B0_prime=$4
    local phase=$5
    
    # Calculate dP/dV and bulk modulus using awk
    awk -v V=$V -v V0=$V0 -v B0=$B0 -v B0_prime=$B0_prime -v PHASE=$phase '
    BEGIN {
        # Calculate dP/dV = (V0^(5/3) * B0 * (5 * V^(2/3) - 7 * V0^(2/3))) / (2 * V^(10/3))
        dPdV = (V0^(5.0/3.0) * B0 * (5.0 * V^(2.0/3.0) - 7.0 * V0^(2.0/3.0))) / (2.0 * V^(10.0/3.0))
        
        # Bulk modulus B(V) = -V * dP/dV
        bulk_modulus = -V * dPdV
        
        # Convert to GPa (assuming input units are consistent)
        bulk_modulus_GPa = bulk_modulus 
        
        printf "Phase: %s\n", PHASE
        printf "Volume: %.6f Å³\n", V
        printf "dP/dV: %.8f\n", dPdV
        printf "Bulk Modulus: %.6f GPa\n", bulk_modulus_GPa
        printf "==========================================\n"
    }'
}

# Main execution
echo "Step 1: Fitting Birch-Murnaghan EOS for both phases..."

# Fit solid phase (column 2)
echo "Fitting solid phase..."
create_gnuplot_fit_script "solid" 2 "solid"
gnuplot fit_solid.gp 2>/dev/null

# Fit liquid phase (column 3)
echo "Fitting liquid phase..."
create_gnuplot_fit_script "liquid" 3 "liquid"
gnuplot fit_liquid.gp 2>/dev/null

echo "Fitting completed!"
echo ""

# Extract fitted parameters
echo "Step 2: Extracting fitted parameters..."

# Read solid parameters
if [[ -f "solid_params.dat" ]]; then
    SOLID_V0=$(grep "V0 =" solid_params.dat | awk '{print $3}')
    SOLID_B0=$(grep "B0 =" solid_params.dat | awk '{print $3}')
    SOLID_B0_PRIME=$(grep "B0_prime =" solid_params.dat | awk '{print $3}')
    
    echo "Solid phase fitted parameters:"
    echo "  V0 = $SOLID_V0"
    echo "  B0 = $SOLID_B0"
    echo "  B0_prime = $SOLID_B0_PRIME"
else
    echo "Error: Could not find solid phase parameters!"
    exit 1
fi

# Read liquid parameters
if [[ -f "liquid_params.dat" ]]; then
    LIQUID_V0=$(grep "V0 =" liquid_params.dat | awk '{print $3}')
    LIQUID_B0=$(grep "B0 =" liquid_params.dat | awk '{print $3}')
    LIQUID_B0_PRIME=$(grep "B0_prime =" liquid_params.dat | awk '{print $3}')
    
    echo "Liquid phase fitted parameters:"
    echo "  V0 = $LIQUID_V0"
    echo "  B0 = $LIQUID_B0"
    echo "  B0_prime = $LIQUID_B0_PRIME"
else
    echo "Error: Could not find liquid phase parameters!"
    exit 1
fi

echo ""

# Step 3: Find volumes for target pressure
echo "Step 3: Finding volumes for target pressure $TARGET_PRESSURE GPa..."

SOLID_VOLUME=$(find_volume_for_pressure $TARGET_PRESSURE $SOLID_V0 $SOLID_B0 $SOLID_B0_PRIME "solid")
LIQUID_VOLUME=$(find_volume_for_pressure $TARGET_PRESSURE $LIQUID_V0 $LIQUID_B0 $LIQUID_B0_PRIME "liquid")

echo "Volume at P=$TARGET_PRESSURE GPa:"
echo "  Solid: $SOLID_VOLUME Å³"
echo "  Liquid: $LIQUID_VOLUME Å³"
echo ""

# Step 4: Calculate bulk modulus
echo "Step 4: Calculating bulk modulus at target pressure..."
echo ""

echo "SOLID PHASE RESULTS:"
calculate_bulk_modulus $SOLID_VOLUME $SOLID_V0 $SOLID_B0 $SOLID_B0_PRIME "Solid"

echo "LIQUID PHASE RESULTS:"
calculate_bulk_modulus $LIQUID_VOLUME $LIQUID_V0 $LIQUID_B0 $LIQUID_B0_PRIME "Liquid"

# Save summary results
cat > summary_results.txt << EOF
Birch-Murnaghan EOS Fitting Results
===================================
Target Pressure: $TARGET_PRESSURE GPa

SOLID PHASE:
- Fitted Parameters:
  V0 = $SOLID_V0
  B0 = $SOLID_B0
  B0_prime = $SOLID_B0_PRIME
- Volume at P=$TARGET_PRESSURE GPa: $SOLID_VOLUME Å³

LIQUID PHASE:
- Fitted Parameters:
  V0 = $LIQUID_V0
  B0 = $LIQUID_B0
  B0_prime = $LIQUID_B0_PRIME
- Volume at P=$TARGET_PRESSURE GPa: $LIQUID_VOLUME Å³

Note: Bulk modulus values are displayed above.
EOF

echo "Results saved to summary_results.txt"
echo "Fit plots saved as solid_fit.png and liquid_fit.png"

# Cleanup temporary files
rm -f fit_solid.gp fit_liquid.gp find_volume.awk
rm -f solid_fitted.dat liquid_fitted.dat

echo ""
echo "Analysis completed successfully!"