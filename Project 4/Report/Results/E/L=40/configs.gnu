source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=40/values(temp)'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=40/accepted configs.png'

set terminal png
set output image

set title "Accepted configurations of a 40x40 lattice"

set xlabel "Temperature"
set ylabel "Accepted configurations

plot source using 1:7 w l title 'Accepted configurations'