set xrange [*:*]
set yrange [*:*]
set title 'F=2.45 m'
set key top left box
set xlabel 'Shell diameter (mm)'
set ylabel 'Integrated Weight (kg)'

plot 'Wtrend.dat' u 2:3 title 'Weight for light (t>0.1mm) shells' w lp,\
'Wtrend.dat' u 2:4 title 'Weight for heavy (t>0.2mm) shells' w lp

set term postscript color 
set out 'thickness_trend2.ps'
replot

set term win
set out

