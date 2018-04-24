set term pngcairo size 1800,1800 color font 'Times Roman, 24'
set output '/home/calat/CLionProjects/DustyShock/cells0.005_CFL0.001_K50000.png'
set style line 1 lt 1 lw 6 lc rgb 'magenta' pt -1
set style line 2 lt 1 lw 3 lc rgb 'black' pt 3 ps 1
set style line 3 lt 1 lw 3 lc rgb 'red' pt 3 ps 1
set style line 4 lt 1 lw 4 lc rgb 'blue' pt 3 ps 1
set style line 5 lt 1 lw 2 lc rgb 'green' pt 3 ps 1
set style line 6 lt 1 lw 2 lc rgb 'black' pt 3 ps 1
set macros
# MACROS
TMARGIN = "set tmargin at screen 0.92; set bmargin at screen 0.66"
MMARGIN = "set tmargin at screen 0.64; set bmargin at screen 0.38"
BMARGIN = "set tmargin at screen 0.36; set bmargin at screen 0.1"
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.49"
RMARGIN = "set lmargin at screen 0.51; set rmargin at screen 0.85"
set multiplot layout 2,2
set key left top
set border ls 6
#set x2label 'Distance'
set x2tics
set ylabel 'Dust density'
set mytics 2 
#set label 1 'explicit, tau=0.0001' at graph 0.4,0.1
set xrange [-1:1]
set x2range [-1:1]
set mx2tics 
#delta = "0.5"
#1
set x2tics mirror
unset xtics
set yrange [0:1.4]
@TMARGIN; @LMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_cellsShock_dust_T0.2_h0.01_tau1e-05_alfa1_beta2_N2250_nu0.1_K50000.dat' using 1:2 title "cells" w l ls 3, "/home/calat/documents/ИК/DustyShock/DustyShockSplashD2G=1T=0.2/dens.dat" using 1:2 title "analytic" w l ls 2

set key right top
#set label 2 'explicit' at graph 0.5,0.1
set y2label 'Dust velocity'
set y2tics
set my2tics 2 
unset ylabel
unset ytics
set y2tics mirror
set yrange [-0.1:1.3]
@TMARGIN; @RMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_cellsShock_dust_T0.2_h0.01_tau1e-05_alfa1_beta2_N2250_nu0.1_K50000.dat' using 1:3 title "cells" w l ls 3, "/home/calat/documents/ИК/DustyShock/DustyShockSplashD2G=1T=0.2/vel.dat" using 1:2 title "analytic" w l ls 2

unset label 1
unset label 2
unset x2label
unset x2tics
#set label 1 'smooth, CFL=0.1' at graph 0.4,0.1
set key left top
set mxtics 2
#set xtics
unset y2label 
unset y2tics
set mytics 2 
set mxtics 1
set ytics
set ylabel 'Gas density'

set format x ""
set xtics
set mxtics
set yrange [0:1.6]
@MMARGIN; @LMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_cellsShock_gas_T0.2_h0.01_tau1e-05_alfa1_beta2_N2250_nu0.1_K50000.dat' using 1:2 title "cells" w l ls 3, "/home/calat/documents/ИК/DustyShock/DustyShockSplashD2G=1T=0.2/dens.dat" using 1:2 title "analytic" w l ls 2

set y2label 'Gas velocity'
set key left top
#set label 2 'smooth' at graph 0.5,0.1
set y2tics

unset ylabel
unset ytics
set y2tics mirror
set yrange [-0.1:1.3]
#set key left bottom
@MMARGIN; @RMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_cellsShock_gas_T0.2_h0.01_tau1e-05_alfa1_beta2_N2250_nu0.1_K50000.dat' using 1:3 title "cells" w l ls 3, "/home/calat/documents/ИК/DustyShock/DustyShockSplashD2G=1T=0.2/vel.dat" using 1:2 title "analytic" w l ls 2

unset format x
unset label 1
unset y2tics
unset y2label
#unset label 2
set xtics
set ylabel 'Energy'
set xlabel 'Distance'
#set key left top
#set label 1 'near, CFL=0.1' at graph 0.4,0.1
set ytics
#unset ylabel
set yrange [1.5:4]
@BMARGIN; @LMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_cellsShock_gas_T0.2_h0.01_tau1e-05_alfa1_beta2_N2250_nu0.1_K50000.dat' using 1:4 title "cells" w l ls 3, "/home/calat/documents/ИК/DustyShock/DustyShockSplashD2G=1T=0.2/ener.dat" using 1:2 title "analytic" w l ls 2


unset ytics
unset ylabel
set y2label 'Pressure'
set y2tics mirror
set yrange [0:1.6]
@BMARGIN; @RMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_cellsShock_gas_T0.2_h0.01_tau1e-05_alfa1_beta2_N2250_nu0.1_K50000.dat' using 1:5 title "cells" w l ls 3, "/home/calat/documents/ИК/DustyShock/DustyShockSplashD2G=1T=0.2/pres.dat" using 1:2 title "analytic" w l ls 2

unset multiplot
