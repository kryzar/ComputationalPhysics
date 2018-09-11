source1 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 1'

source2 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 2'

source3 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 3'

source4 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 4'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/énergies positive.png'

set terminal png enhanced
set output image

set title "Mean energy of a 2x2 positive lattice"
set xrange[0:1000000]
set xlabel "Monte Carlo cycles"
set ylabel "<E>"

plot source1 using 1:2 w l title 'Mean energy 1', source2 using 1:2 w l title 'Mean energy 2', source3 using 1:2 w l title 'Mean energy 3', source4 using 1:2 w l title 'Mean energy 4'
