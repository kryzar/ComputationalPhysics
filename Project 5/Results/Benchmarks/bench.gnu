reset

source = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Benchmarks/onedim.txt"

image = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 5/Calypso/Results/Benchmarks/bench dx**2"

set terminal png
set output image
## set xrange[0.01:0.000001]

set logscale xy
set xlabel "dx^2"
set ylabel "log(|relative error|)"

plot source using 2:9 w lp title "Crank-Nicolson scheme"
