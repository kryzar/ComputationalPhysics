source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=40/values(temp)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=40/khi.png'

set terminal png
set output image

set title "Susceptibility of a 40x40 lattice"

set xlabel "Temperature"
set ylabel "Suceptibility"

plot source using 1:5 w l title 'Suceptibility'