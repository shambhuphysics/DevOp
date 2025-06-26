#!/bin/bash

# Birch-Murnaghan EOS Fitting - Dynamic Parameters Version
INPUT_FILE="VP.dat"
TARGET_PRESSURE=0.5

# Parse arguments
while getopts "f:p:h" opt; do
    case $opt in
        f) INPUT_FILE="$OPTARG" ;;
        p) TARGET_PRESSURE="$OPTARG" ;;
        h) echo "Usage: $0 [-f FILE] [-p PRESSURE]"; exit 0 ;;
    esac
done

# Check if file exists
[[ ! -f "$INPUT_FILE" ]] && { echo "Error: $INPUT_FILE not found!"; exit 1; }

echo "Birch-Murnaghan EOS Fitting - Target: $TARGET_PRESSURE GPa"

# Function to estimate initial parameters from data
estimate_initial_params() {
    local column=$1
    awk -v col=$column '
    !/^#/ && NF >= 3 {
        V[++n] = $1; P[n] = $col
        V_sum += $1; P_sum += $col
        if (n == 1 || $1 > V_max) V_max = $1
        if (n == 1 || $1 < V_min) V_min = $1
        if (n == 1 || $col > P_max) P_max = $col
        if (n == 1 || $col < P_min) P_min = $col
    }
    END {
        if (n < 3) { print "11000 100 4.5"; exit }
        
        # Estimate V0 as extrapolated volume at P=0
        V0_est = V_max + (V_max - V_min) * 0.2
        
        # Estimate B0 from pressure range and volume compression
        V_range = V_max - V_min
        P_range = P_max - P_min
        if (V_range > 0 && P_range > 0) {
            # Rough estimate: B0 ~ -V * dP/dV
            dPdV_est = -P_range / V_range
            B0_est = (V_max + V_min) / 2 * dPdV_est
            if (B0_est < 10) B0_est = 50
            if (B0_est > 500) B0_est = 200
        } else {
            B0_est = 100
        }
        
        # B0_prime typically ranges 3-6 for most materials
        B0_prime_est = 4.5
        
        printf "%.0f %.0f %.1f\n", V0_est, B0_est, B0_prime_est
    }' "$INPUT_FILE"
}

# Function to get data range for bisection bounds
get_data_range() {
    awk '!/^#/ && NF >= 3 {
        if (NR == 1 || $1 < min_V) min_V = $1
        if (NR == 1 || $1 > max_V) max_V = $1
    }
    END {
        # Extend range by 20% on each side for safety
        range = max_V - min_V
        extension = range * 0.2
        printf "%.0f %.0f\n", min_V - extension, max_V + extension
    }' "$INPUT_FILE"
}

# Function to fit and analyze phase
analyze_phase() {
    local phase=$1 column=$2
    
    # Get initial parameter estimates
    read V0_init B0_init B0_prime_init <<< "$(estimate_initial_params $column)"
    
    # Get volume range for bisection
    read V_low_range V_high_range <<< "$(get_data_range)"
    
    echo "Initial estimates for ${phase}: V0=$V0_init, B0=$B0_init, B0'=$B0_prime_init"
    
    # Create gnuplot fit script with dynamic initial values
    cat > fit_${phase}.gp << EOF
        birch(V, V0, B0, B0_prime) = (3.0/2.0) * B0 * ((V0/V)**(7.0/3.0) - (V0/V)**(5.0/3.0)) * (1 + (3.0/4.0) * (B0_prime - 4) * ((V0/V)**(2.0/3.0) - 1))
        V0 = $V0_init; B0 = $B0_init; B0_prime = $B0_prime_init
        fit birch(x, V0, B0, B0_prime) '$INPUT_FILE' using 1:$column via V0, B0, B0_prime
        set print '${phase}_params.dat'
        print "V0 = ", V0
        print "B0 = ", B0
        print "B0_prime = ", B0_prime
        set print
EOF
    
    # Run gnuplot fit
    gnuplot fit_${phase}.gp 2>/dev/null
    
    # Extract parameters
    if [[ -f "${phase}_params.dat" ]]; then
        V0=$(grep "V0 =" ${phase}_params.dat | awk '{print $3}')
        B0=$(grep "B0 =" ${phase}_params.dat | awk '{print $3}')
        B0_prime=$(grep "B0_prime =" ${phase}_params.dat | awk '{print $3}')
        
        # Find volume using bisection method with dynamic bounds
        VOLUME=$(awk -v TARGET_P=$TARGET_PRESSURE -v V0_VAL=$V0 -v B0_VAL=$B0 -v B0_PRIME_VAL=$B0_prime \
                     -v V_LOW=$V_low_range -v V_HIGH=$V_high_range '
        function birch_pressure(V, V0, B0, B0_prime) {
            return (3.0/2.0) * B0 * ((V0/V)^(7.0/3.0) - (V0/V)^(5.0/3.0)) * (1 + (3.0/4.0) * (B0_prime - 4) * ((V0/V)^(2.0/3.0) - 1))
        }
        function abs(x) { return x < 0 ? -x : x }
        BEGIN {
            target = TARGET_P; V0 = V0_VAL; B0 = B0_VAL; B0_prime = B0_PRIME_VAL
            V_low = V_LOW; V_high = V_HIGH
            tolerance = 1e-6; max_iter = 1000
            
            # Ensure bounds make sense
            if (V_low <= 0) V_low = V0 * 0.5
            if (V_high <= V_low) V_high = V0 * 1.5
            
            for (i = 1; i <= max_iter; i++) {
                V_mid = (V_low + V_high) / 2.0
                P_mid = birch_pressure(V_mid, V0, B0, B0_prime)
                
                if (abs(P_mid - target) < tolerance) break
                
                P_low = birch_pressure(V_low, V0, B0, B0_prime)
                if ((P_low - target) * (P_mid - target) < 0) {
                    V_high = V_mid
                } else {
                    V_low = V_mid
                }
            }
            
            # Calculate bulk modulus
            dPdV = (V0^(5.0/3.0) * B0 * (5.0 * V_mid^(2.0/3.0) - 7.0 * V0^(2.0/3.0))) / (2.0 * V_mid^(10.0/3.0))
            bulk_modulus = -V_mid * dPdV
            
            printf "%.2f %.2f\n", V_mid, bulk_modulus
        }')
        
        read vol bulk <<< "$VOLUME"
        
        echo "${phase^} Phase:"
        echo "  V0=$V0, B0=$B0, B0'=$B0_prime"
        echo "  Volume at $TARGET_PRESSURE GPa: $vol Å³"
        echo "  Bulk Modulus: $bulk GPa"
        echo
    else
        echo "Error: Could not fit ${phase} phase!"
    fi
    
    # Cleanup
    rm -f fit_${phase}.gp ${phase}_params.dat
}

# Analyze both phases
analyze_phase "solid" 2
analyze_phase "liquid" 3

echo "Analysis completed!"
