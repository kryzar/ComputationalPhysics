path_source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Program/Results/arma v jacobi.txt'

path_output = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Report/Plots/Times/plot.tex'

set xlabel "n"
set ylabel "times (s)"

set terminal latex
set output path_output


plot path_source using 1:2 w l title 'Jacobi', path_source using 1:3 w l title 'Armadillo'