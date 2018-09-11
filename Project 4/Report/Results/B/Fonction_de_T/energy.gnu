source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Fonction de T/simulation'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Fonction de T/energy.png'

set terminal png
set output image

set title "Mean energy of a 2x2 lattice"

set xlabel "Temperature"
set ylabel "<E>"

plot source using 1:2 w l title 'Mean energy'