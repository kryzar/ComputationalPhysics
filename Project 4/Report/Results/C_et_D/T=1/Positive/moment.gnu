source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C et D/T=1/Positive/values(cycle)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C et D/T=1/Positive/moment.png'

set terminal png
set output image

set title "Mean magnetic moment of a 20x20 positive lattice, T=1"
set xlabel "Monte Carlo cycles"
set ylabel "<|M|>"

plot source using 1:3 w l title 'Mean magnetic moment'