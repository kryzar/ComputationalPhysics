source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/values(cycle)'

set xlabel "Monte Carlo cycles"
set ylabel "Accepted configurations"

plot source using 1:7 w l title 'Accepted configurations'