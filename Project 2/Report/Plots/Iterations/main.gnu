path_source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Program/Results/arma v jacobi.txt'

path_output = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Report/Plots/Iterations/plot.tex'

set xlabel "n"

set terminal latex
set output path_output


plot path_source using 1:5 w l title 'iterations'