plot "Verify/Output_vel_Re400_Eps-4.csv" u 3:($2/128)  every 129::66 with linespoint
replot "Verify/Ghia_results/ghia-400u.dat" u 1:2
set xlabel 'u'
set ylabel 'y'
pause -1
