# Advanced Molecular Dynamics Live Dashboard
# Usage: gnuplot dashboard.gp
# Data sources: OSZICAR, density.dat, XDATCAR (optional)

reset
set term x11 size 1200,900 enhanced
set multiplot layout 3,3

# Global settings
set grid alpha 0.3
set key outside right
set tics font ",10"
set label font ",12"

# Header with system info
set label 1 at screen 0.5,0.97 center "{/Bold=16 MD Simulation Dashboard}" 
set label 2 at screen 0.5,0.94 center "{/Bold=12 ".system("pwd")."}"
set label 3 at screen 0.02,0.91 "{/Bold Timestamp: ".strftime("%Y-%m-%d %H:%M:%S")."}" font ",10"

# Plot 1: Temperature vs Time
set title '{/Bold=14 Temperature Evolution}' offset 0,-0.5
set xlabel '{/Bold Step}' font ",12"
set ylabel '{/Bold Temperature (K)}' font ",12"
set format y "%.0f"
plot 'OSZICAR' u 1:3 w l lw 2 lc rgb "#FF4444" title "T_{inst}", \
     '' u 1:3 smooth bezier lw 1 lc rgb "#AA0000" title "T_{smooth}"

# Plot 2: Energy Components  
set title '{/Bold=14 Energy Components}' offset 0,-0.5
set xlabel '{/Bold Step}' font ",12"
set ylabel '{/Bold Energy (eV)}' font ",12"
set format y "%.2f"
plot 'OSZICAR' u 1:7 w l lw 2 lc rgb "#00AA00" title "Total E", \
     '' u 1:4 w l lw 2 lc rgb "#0066CC" title "Kinetic E" if exists

# Plot 3: Constant of Motion (Energy Conservation)
set title '{/Bold=14 Energy Conservation}' offset 0,-0.5
set xlabel '{/Bold Step}' font ",12" 
set ylabel '{/Bold Conserved Quantity}' font ",12"
set format y "%.4f"
stats 'OSZICAR' u 5 nooutput
plot 'OSZICAR' u 1:5 w l lw 2 lc rgb "#4444FF" title "Const. Motion", \
     STATS_mean w l lw 1 dt 2 lc rgb "#000000" title sprintf("Mean: %.4f", STATS_mean)

# Plot 4: Density Profile
set title '{/Bold=14 Density Distribution}' offset 0,-0.5
set xlabel '{/Bold Position (slabs)}' font ",12"
set ylabel '{/Bold Atom Count}' font ",12"
set yrange [0:*]
set format y "%.0f"
plot 'density.dat' u 1:2 w filledcurves x1 lc rgb "#AA44AA" fillstyle solid 0.3 title "Density", \
     '' u 1:2 w l lw 2 lc rgb "#660066" notitle

# Plot 5: Pressure (if available)
set title '{/Bold=14 Pressure Evolution}' offset 0,-0.5
set xlabel '{/Bold Step}' font ",12"
set ylabel '{/Bold Pressure (GPa)}' font ",12"
if (system("test -f OSZICAR && awk 'NF>=9{exit 0} END{exit 1}' OSZICAR") == 0) {
    plot 'OSZICAR' u 1:8 w l lw 2 lc rgb "#FF8800" title "Pressure"
} else {
    set label 10 at graph 0.5,0.5 center "Pressure data\nnot available" font ",14"
    plot NaN notitle
    unset label 10
}

# Plot 6: Volume/Cell Evolution  
set title '{/Bold=14 Cell Volume}' offset 0,-0.5
set xlabel '{/Bold Step}' font ",12"
set ylabel '{/Bold Volume (Å³)}' font ",12"
if (system("test -f OSZICAR && awk 'NF>=10{exit 0} END{exit 1}' OSZICAR") == 0) {
    plot 'OSZICAR' u 1:9 w l lw 2 lc rgb "#00CCCC" title "Volume"
} else {
    set label 11 at graph 0.5,0.5 center "Volume data\nnot available" font ",14"  
    plot NaN notitle
    unset label 11
}

# Plot 7: Temperature Distribution (Histogram)
set title '{/Bold=14 T Distribution}' offset 0,-0.5
set xlabel '{/Bold Temperature (K)}' font ",12"
set ylabel '{/Bold Frequency}' font ",12"
set boxwidth 10
bin_width = 10
bin(x) = bin_width * floor(x/bin_width) + bin_width/2.0
plot 'OSZICAR' u (bin($3)):(1.0) smooth freq w boxes lc rgb "#FF6666" fillstyle solid 0.5 title "T Hist"

# Plot 8: Phase Space (T vs E)
set title '{/Bold=14 Phase Space (T-E)}' offset 0,-0.5  
set xlabel '{/Bold Temperature (K)}' font ",12"
set ylabel '{/Bold Energy (eV)}' font ",12"
plot 'OSZICAR' u 3:7 w p pt 7 ps 0.3 lc rgb "#888888" title "T-E trajectory"

# Plot 9: System Status & Statistics
set title '{/Bold=14 System Status}' offset 0,-0.5
unset xlabel
unset ylabel
unset xtics
unset ytics
unset border

# Calculate and display statistics
stats 'OSZICAR' u 3 name "T" nooutput
stats 'OSZICAR' u 7 name "E" nooutput  
stats 'OSZICAR' u 5 name "C" nooutput

set label 20 at graph 0.05,0.9 "{/Bold=12 Current Statistics:}" left
set label 21 at graph 0.05,0.8 sprintf("T: %.1f ± %.1f K", T_mean, T_stddev) left
set label 22 at graph 0.05,0.7 sprintf("E: %.3f ± %.3f eV", E_mean, E_stddev) left  
set label 23 at graph 0.05,0.6 sprintf("Drift: %.2e", C_stddev/abs(C_mean)*100) left
set label 24 at graph 0.05,0.5 sprintf("Steps: %.0f", T_records) left
set label 25 at graph 0.05,0.4 sprintf("Runtime: %.1f min", T_records*0.001) left

# System info
set label 26 at graph 0.05,0.25 "{/Bold=12 Files:}" left
set label 27 at graph 0.05,0.15 system("ls -la OSZICAR CONTCAR 2>/dev/null | tail -2 | awk '{print $9\": \"$6\" \"$7\" \"$8}'") left

plot NaN notitle

# Cleanup labels
unset label 20; unset label 21; unset label 22; unset label 23
unset label 24; unset label 25; unset label 26; unset label 27

unset multiplot

# Auto-refresh settings
pause 2
reread
