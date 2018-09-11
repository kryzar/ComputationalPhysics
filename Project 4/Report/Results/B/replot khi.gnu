source5 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 1'

source6 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 2'

source7 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 3'

source8 = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/T=3/Lattice positive 4'

image = '/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/B/Après correction/khi positive.png'

set terminal png enhanced
set output image

set title "Susceptibility of a 2x2 positive lattice"
set xlabel "Monte Carlo cycles"
set ylabel "Khi"

plot source5 using 1:6 w l title 'Susceptibility 1', source6 using 1:6 w l title 'Susceptibility 2', source7 using 1:6 w l title 'Susceptibility 3', source8 using 1:6 w l title 'Susceptibility 4'