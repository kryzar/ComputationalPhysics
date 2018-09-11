source1 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 1'

source2 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 2'

source3 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 3'

source4 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 4'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/énergies random.png'

set terminal png enhanced
set output image

set xlabel "Monte Carlo cycles"
set ylabel "<E>"

plot source1 using 1:2 w l title 'Mean energy 1', source2 using 1:2 w l title 'Mean energy 2', source3 using 1:2 w l title 'Mean energy 3', source4 using 1:2 w l title 'Mean energy 4'
