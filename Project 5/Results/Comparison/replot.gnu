reset

set size ratio -1
imp = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Tests/imp"

exp = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Tests/exp"

cn = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Tests/cn"

ana = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Tests/ana"

set terminal png
set output "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Tests/replot.png"

set title "Comparison of the algorithms, t=0.200000"
set xlabel "x"
set ylabel "u(x,t)"

plot exp using 1:2 w l title "Explicit", imp using 1:2 w l title "Implicit", cn using 1:2 w l title "Crank-Nicolson", ana using 1:2 w l title "Analytical solution"
