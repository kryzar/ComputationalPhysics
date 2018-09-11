reset

set size ratio -1
data = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/One dimension/t = 0,00001/dx = 0,1/CN scheme/results"

set terminal png
set output "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/One dimension/t = 0,00001/dx = 0,1/CN scheme/map 0.000010.png"
set title "Diffusion equation in one dimension, Crank-Nicolson scheme, t=0.000010"
set xlabel "x"
set ylabel "u(x,t)"

plot "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/One dimension/t = 0,00001/dx = 0,1/CN scheme/results" using 1:2 w l notitle
