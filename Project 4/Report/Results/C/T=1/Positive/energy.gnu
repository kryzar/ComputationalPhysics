source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=1/Positive/values(cycle)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=1/Positive/energy.png'

set terminal png
set output image

set title "Mean energy of a 20x20 positive lattice, T=1"
set xlabel "Monte Carlo cycles"
set ylabel "<E>"

plot source using 1:2 w l title 'Mean energy'