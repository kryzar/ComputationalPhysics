reset

set size ratio -1
set title "Diffusion equation in one dimension, implicit scheme, t=0.300000"
set xlabel "x"
set ylabel "u(x,t)"

plot "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Benchmarks/results" using 1:2 w l notitle
