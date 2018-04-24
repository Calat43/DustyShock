set term pngcairo size 1800,1200 color font 'Times Roman, 24'
set output '/home/calat/CLionProjects/DustyShock/tmp.png'
set style line 1 lt 1 lw 3 lc rgb 'magenta' pt 2 ps 1
set style line 2 lt 1 lw 2 lc rgb 'black' pt -1
set style line 3 lt 1 lw 3 lc rgb 'green' pt 2 ps 1
set style line 4 lt 1 lw 4 lc rgb 'magenta' pt -1
set style line 5 lt 1 lw 2 lc rgb 'black' pt -1
set style line 6 lt 1 lw 3 lc rgb 'blue' pt 2 ps 1
set macros
# MACROS
TMARGIN = "set tmargin at screen 0.9; set bmargin at screen 0.51"
BMARGIN = "set tmargin at screen 0.49; set bmargin at screen 0.1"
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.49"
RMARGIN = "set lmargin at screen 0.51; set rmargin at screen 0.85"
set multiplot layout 2,2
set key right top
#set key at -0.27, 0.0001
set border ls 5
set x2label 'Distance'
set x2tics
#set ylabel 'Dust velocity'
set mytics 2 
#set label 1 't_{stop}=0.002, {/Symbol e}=1' at graph 0.4, 0.2
#set label 3 'h=0.05, {/Symbol t}=0.005' at graph 0.4, 0.1
#set yrange [-0.000105:0.000105]
set xrange [0:1]
set x2range [0:1]
set mx2tics 1
set x2tics mirror
unset xtics
@TMARGIN; @LMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.01_tau0.001_alfa1_beta2_N2250_nu0.1.dat' u 1:2 axes x2y1 title "rho, N=2250" w l ls 1, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.005_tau0.0005_alfa1_beta2_N4500_nu0.1.dat' u 1:2 axes x2y1 title "N=4500" w l ls 6, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.0025_tau0.00025_alfa1_beta2_N9000_nu0.1.dat' u 1:2 axes x2y1 title "N=9000" w l ls 3, '/home/calat/documents/ИК/DustyShock/analyticT=0.2/dens.dat' u 1:2 axes x2y1 title "analytic" w l ls 2 
unset label 1
#set key at 1, 0.0001


unset label 3
set key left top
#set label 2 't_{stop}=0.002, {/Symbol e}=1' at graph 0.4, 0.2
#set label 4 'h=0.05, {/Symbol t}=0.001' at graph 0.4, 0.1
#set y2label 'Dust velocity'
set y2tics
set y2range [-0.52:3.7]
set my2tics 2 
unset ylabel
unset ytics
set y2tics mirror
@TMARGIN; @RMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.01_tau0.001_alfa1_beta2_N2250_nu0.1.dat' u 1:3 axes x2y2 title "vel, N=2250" w l ls 1, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.005_tau0.0005_alfa1_beta2_N4500_nu0.1.dat' u 1:3 axes x2y2 title "N=4500" w l ls 6, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.0025_tau0.00025_alfa1_beta2_N9000_nu0.1.dat' u 1:3 axes x2y2 title "N=9000" w l ls 3, '/home/calat/documents/ИК/DustyShock/analyticT=0.2/vel.dat' u 1:2 axes x2y2 title "analytic" w l ls 2 
unset label 2
unset label 4
unset x2label
unset x2tics
set key left top
#set key at 0, 1.00015

set mxtics 2
set xtics
unset y2label 
unset y2tics
#set yrange [0.99985:1.00015]
#set y2range [0.99985:1.00015]
set mytics 2 
set mxtics 1
set ytics
#set ylabel 'Gas density'
set xlabel 'Distance'
@BMARGIN; @LMARGIN
plot '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.01_tau0.001_alfa1_beta2_N2250_nu0.1.dat' u 1:4 axes x2y1 title "ener, N=2250" w l ls 1, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.005_tau0.0005_alfa1_beta2_N4500_nu0.1.dat' u 1:4 axes x2y1 title "N=4500" w l ls 6, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.0025_tau0.00025_alfa1_beta2_N9000_nu0.1.dat' u 1:4 axes x2y1 title "N=9000" w l ls 3, '/home/calat/documents/ИК/DustyShock/analyticT=0.2/ener.dat' u 1:2 axes x2y1 title "analytic" w l ls 2 
#set key at -0.27, 1.00015
#set y2label 'Gas density'
set y2tics
set y2range [0.07:1.11]
unset ylabel
unset ytics
set y2tics mirror
@BMARGIN; @RMARGIN
#set key left bottom
plot '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.01_tau0.001_alfa1_beta2_N2250_nu0.1.dat' u 1:5 axes x2y2 title "pres, N=2250" w l ls 1, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.005_tau0.0005_alfa1_beta2_N4500_nu0.1.dat' u 1:5 axes x2y2 title "N=4500" w l ls 6, '/home/calat/CLionProjects/DustyShock/im_shock_T0.2_h0.0025_tau0.00025_alfa1_beta2_N9000_nu0.1.dat' u 1:5 axes x2y2 title "N=9000" w l ls 3, '/home/calat/documents/ИК/DustyShock/analyticT=0.2/pres.dat' u 1:2 axes x2y2 title "analytic" w l ls 2 
unset multiplot
