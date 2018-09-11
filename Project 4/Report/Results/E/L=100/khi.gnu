source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=100/values(temp)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=100/khi.png'

set terminal png
set output image

set title "Susceptibility of a 20x20 lattice"

set xlabel "Temperature"
set ylabel "Susceptibility"

plot source using 1:6 w l title 'Susceptibility'