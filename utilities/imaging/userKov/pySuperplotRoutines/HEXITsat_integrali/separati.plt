set ylabel 'Effective Area (cm^2) '
set xlabel 'Shell Number from outer'
set title 'Integrated Effective Area (cm^2) @ 10 keV for F10D350ff010_thsx, F=, ff=, Dmax= configuration'
plot 'area10_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2)  x 1.0' w lp
pause-1

set ylabel 'Effective Area (cm^2) '
set xlabel 'Shell Number from outer'
set title 'Integrated Effective Area (cm^2) @ 30 keV for F10D350ff010_thsx, F=, ff=, Dmax= configuration'
plot 'area30_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2)  x 1.0' w lp
pause-1

set ylabel 'Effective Area (cm^2) '
set xlabel 'Shell Number from outer'
set title 'Integrated Effective Area (cm^2) @ 60 keV for F10D350ff010_thsx, F=, ff=, Dmax= configuration'
plot 'area60_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2)  x 1.0' w lp
pause-1

set ylabel 'Effective Area (cm^2) '
set xlabel 'Shell Number from outer'
set title 'Integrated Effective Area (cm^2) @ 70 keV for F10D350ff010_thsx, F=, ff=, Dmax= configuration'
plot 'area70_F10D350ff010_thsx.txt' u 1:3 title 'Effective Area (cm^2)  x 1.0' w lp
pause-1

set ylabel 'Colleting Area (cm^2) '
set xlabel 'Shell Number from outer'
set title 'Integrated Colleting Area (cm^2)  for F10D350ff010_thsx, F=, ff=, Dmax= configuration'
plot 'Acoll_F10D350ff010_thsx.txt' u 1:3 title 'Colleting Area (cm^2)  x 1.0' w lp
pause-1

set ylabel 'Weight (kg)'
set xlabel 'Shell Number from outer'
set title 'Integrated Weight (kg) for F10D350ff010_thsx, F=, ff=, Dmax= configuration'
plot 'Weight_F10D350ff010_thsx.txt' u 1:3 title 'Weight (kg) x 1.0' w lp
pause-1

