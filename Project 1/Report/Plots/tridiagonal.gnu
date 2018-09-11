path = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 1/Alcyonide/Programs/Tridiagonal case/Data files/results n=1000'

set xlabel "xi"

set terminal latex
set output '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 1/Alcyonide/Report/plots/tri1000.tex'


plot path using 2:3 w l title 'numerical solution', path using 2:4 w l title 'exact solution'