density_profile() {
    local ndim=400
    local infile="${1:-CONTCAR}"
    local outfile="${2:-density.dat}"
    
    awk -v ndim="$ndim" '
    BEGIN {
        dx = 1.0/ndim
        for(j=1; j<=ndim; j++) d[j] = 0
    }
    NR==2 { volume = $1 }
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
density_profile  
