source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=100/values(temp)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=100/energy.png'

set terminal png
set output image

set title "Mean energy of a 20x20 lattice"

set xlabel "Temperature"
set ylabel "<E>"

plot source using 1:2 w l title 'Mean energy'