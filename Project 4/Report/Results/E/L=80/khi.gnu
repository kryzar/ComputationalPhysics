source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=80/values(temp)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=80/khi.png'

set terminal png
set output image

set title "Susceptibility of a 80x80 lattice"

set xlabel "Temperature"
set ylabel "Susceptibility"

plot source using 1:5 w l title 'Susceptibility'