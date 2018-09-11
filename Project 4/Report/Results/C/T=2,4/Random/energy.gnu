source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=2,4/Random/values(cycle)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=2,4/Random/energy.png'

set terminal png
set output image

set title "Mean energy of a 20x20 random lattice, T=2,4"
set xlabel "Monte Carlo cycles"
set ylabel "<E>"

plot source using 1:2 w l title 'Mean energy'