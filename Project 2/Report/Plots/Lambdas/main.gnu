path_source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Report/Plots/Lambdas/plot.txt'

path_output = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Report/Plots/Lambdas/plot.tex'

set xlabel '$\omega$'
set ylabel '${\lambda}_{n}$'

set terminal latex
set output path_output

plot path_source using 1:2 w lp title '${\lambda}_{0}$', path_source using 1:3 w lp title '${\lambda}_{1}$', path_source using 1:4 w lp title '${\lambda}_{2}$' 