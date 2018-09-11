path = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 1/Alcyonide/Programs/Tridiagonal case/Data files/errors'

set xlabel "n"

set terminal latex
set output '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 1/Alcyonide/Report/plots/errors.tex'


plot path using 1:2 w l title 'relative error'