set term x11 persist size 1200,800

set title "Energy and Pressure Comparison"
set xlabel "Time step"
set ylabel "Energy (eV)/Pressure (kBr)"
set border linewidth 1.5

set style line 1 lc rgb 'red' lt 1 lw 2
set style line 2 lc rgb 'blue' lt 1 lw 2
set style line 3 lc rgb 'red' lt 1 lw 2 pt 7
set style line 4 lc rgb 'blue' lt 1 lw 2 pt 7

set key samplen 2 spacing 1.25 font ",14"

set multiplot layout 2,1 spacing 0.15, 0.25

# First subplot for energy comparison
set title "Energy Fit" font ",20"
plot 'data_temp' u 1:2 w l ls 1 ti 'Energy_EAM', \
     '' u 1:3 w l ls 2 ti 'Energy_DFT'

# Second subplot for pressure comparison
set title "Pressure Fit" font ",20"
plot 'data_temp' u 1:5 w l ls 3 ti 'Pressure_EAM', \
     '' u 1:6 w l ls 4 ti 'Pressure_DFT'

unset multiplot
pause 0.01

reread

