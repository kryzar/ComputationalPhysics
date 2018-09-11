source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/values(cycle)'

set xlabel "Monte Carlo cycles"
set ylabel "<|M|>"

plot source using 1:3 w l title 'Magnetic moment'