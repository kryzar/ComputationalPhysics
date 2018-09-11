source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Fonction de T/simulation'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Fonction de T/cv.png'

set terminal png
set output image

set title "Mean specific heat of a 2x2 lattice"

set xlabel "Temperature"
set ylabel "Specific heat"

plot source using 1:5 w l title 'Mean specific heat'