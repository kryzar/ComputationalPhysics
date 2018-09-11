source5 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 1'

source6 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 2'

source7 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 3'

source8 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 4'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/cv random.png'

set terminal png enhanced
set output image

set xlabel "Monte Carlo cycles"
set ylabel "Cv"

plot source5 using 1:5 w l title 'Specific heat 1', source6 using 1:5 w l title 'Specific heat 2', source7 using 1:5 w l title 'Specific heat 3', source8 using 1:5 w l title 'Specific heat 8'