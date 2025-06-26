# Function to estimate the Pressure
calculate_average_pressure() {
    local outcar_file="${1:-OUTCAR}"
    
    # Check if file exists
    if [[ ! -f "$outcar_file" ]]; then
        echo "Error: $outcar_file file not found."
        return 1
    fi
    
    # Extract lines containing 'plus kinetic' and calculate averages
    grep 'plus kinetic' "$outcar_file" | awk '{
        # Calculate average of columns 3, 4, 5
        avg = ($3 + $4 + $5) / 3
        # Sum averages for final calculation
        sum += avg
        count++
    } END {
        # Calculate mean of averages, divide by 10, and print result
        if (count > 0) {
            mean = (sum / count) / 10
            printf "Pressure at this stage: %.2f GPa\n", mean
        } else {
            print "No plus kinetic lines found in file"
        }
    }' 
}

calculate_average_pressure