source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=1/Positive/values(cycle)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=1/Positive/configs.png'

set terminal png
set output image

set title "Number of accepted configurations of a 20x20 positive lattice, T=1"
set xlabel "Monte Carlo cycles"
set ylabel "Accepted configurations"

plot source using 1:7 w l title 'Accepted configurations'