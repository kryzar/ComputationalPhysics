reset

set size ratio -1
set title "Diffusion equation in one dimension, explicit scheme, t=0.000010"
set xlabel "x"
set ylabel "u(x,t)"

plot "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/One dimension/dx = 0,1/Explicit scheme/results" using 1:2 w l notitle
