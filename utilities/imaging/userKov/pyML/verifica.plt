plot 'interp.txt' u 1:2 w lp,\
'iridio.dat' u 1:(1-$2) w l
pause -1

plot 'interp.txt' u 1:3 w lp,\
'iridio.dat' u 1:3 w l
pause -1

plot 'IrEner.txt' u 1:2 w lp,\
'kovIr.dat' u 1:2 w p