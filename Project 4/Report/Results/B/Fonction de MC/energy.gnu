source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/values(cycle)'

set xlabel "Monte Carlo cycles"
set ylabel "<E>"

plot source using 1:2 w l title 'Mean energy'