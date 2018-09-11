source5 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 5'

source6 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 6'

source7 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 7'

source8 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/lattice random 8'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/susceptibility random.png'

set terminal png enhanced
set output image

set xlabel "Monte Carlo cycles"
set ylabel "Khi"

plot source5 using 1:6 w l title 'Suceptibility 1', source6 using 1:6 w l title 'Suceptibility 2', source7 using 1:6 w l title 'Suceptibility 3', source8 using 1:6 w l title 'Suceptibility 8'