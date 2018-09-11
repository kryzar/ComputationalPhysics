reset

set size ratio -1
data = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/One dimension/t = 0,2/dx = 0,1/Implicit scheme/results"

set terminal png
set output "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/One dimension/t = 0,2/dx = 0,1/Implicit scheme/map 0.200000.png"
set title "Diffusion equation in one dimension, implicit scheme, t=0.200000"
set xlabel "x"
set ylabel "u(x,t)"

plot "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/One dimension/t = 0,2/dx = 0,1/Implicit scheme/results" using 1:2 w l notitle
