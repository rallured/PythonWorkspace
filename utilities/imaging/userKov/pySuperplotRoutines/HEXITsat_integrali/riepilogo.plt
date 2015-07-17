set title 'F10D350ff010_thsx, F=, ff=, Dmax='
set key top left box
set xlabel 'Shell Number from outer'
set yrange [0:1000] 
set ylabel 'Effective Area (cm^2) '
plot 'area10_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2) @ 10 keV x 1.0' w lp,\
'area30_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2) @ 30 keV x 1.0' w lp,\
'area60_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2) @ 60 keV x 1.0' w lp,\
'area70_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2) @ 70 keV x 1.0' w lp,\
'Acoll_F10D350ff010_thsx.txt' u 1:3 title 'Colleting Area (cm^2)  x 1.0' w lp
pause-1

set terminal postscript color
set output '000_Acoll_F10D350ff010_thsx.ps'
replot
set term win
set out

set terminal png
set output '000_Acoll_F10D350ff010_thsx.png'
replot
set term win
set out

set yrange [0:100] 
set ylabel 'Colleting Area (cm^2) '
plot 'Acoll_F10D350ff010_thsx.txt' u 1:3 title 'Colleting Area (cm^2)  x 1.0' w lp
pause-1

set terminal postscript color
set output '001_Acoll_F10D350ff010_thsx.ps'
replot
set term win
set out

set terminal png
set output '001_Acoll_F10D350ff010_thsx.png'
replot
set term win
set out

