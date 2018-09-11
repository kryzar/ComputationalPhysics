reset

source = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=2,4/Positive/energies proba'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/C/T=2,4/Positive/energies proba.png'

set term png truecolor
set output image

set title "Probability of finding a given energy in the Monte-Carlo cycles"
set xlabel "Energy
set ylabel "Probability"
set grid
set style fill transparent solid 0.5 noborder
plot source u 1:2 w boxes lc rgb "green" notitle
