source5 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 1'

source6 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 2'

source7 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 3'

source8 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 4'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/moments random.png'

set terminal png enhanced
set output image

set title "Mean magnetic moment of a 2x2 positive lattice"
set xlabel "Monte Carlo cycles"
set ylabel "<|M|>"

plot source5 using 1:3 w l title 'Mean magnetization 1', source6 using 1:3 w l title 'Mean magnetization 2', source7 using 1:3 w l title 'Mean magnetization 3', source8 using 1:3 w l title 'Mean magnetization 4'