source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=60/values(temp)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=60/moment.png'

set terminal png
set output image

set title "Mean absolute magnetic moment of a 60x60 lattice"

set xlabel "Temperature"
set ylabel "<|M|>"

plot source using 1:3 w l title 'Mean magnetic moment'