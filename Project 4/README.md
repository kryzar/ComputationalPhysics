# Ising Peasy <br> Project 4 from the course FYS3150 of University of Oslo, Autumn 2017

This program was made to compute expected values such as the mean energy, the mean absolute value of the magnetization, the susceptibility and the specific heat of a two-dimensional lattice in the Ising Model.

[![moments_random.png](https://s17.postimg.org/t27fie83j/moments_random.png)](https://postimg.org/image/5b820a7wb/)

## Tutorial

### Basic computations

This program is object-oriented and is based on the class `lattice`. To declare or initialize a 10x10 squared lattice called `system` with initial temperature 1 (normalized with the Boltzman constant and the J constant), just :

```cpp
#include "lattice.hpp"

int main(int argc, char* argv[])
{
  lattice system(10, 1.); //  10x10 lattice with temperature 1.

  return 0;
}
```

Note that all the spins are 1 by default. After that, you can compute simple operations on the lattice and see the output :

```cpp
#include "lattice.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
    lattice system(10, 1.); //  10x10 lattice with temperature 1.

    system.print();         //  prints the lattice
    system.randomize();     //  randomize all the spins
    system.print();
    return 0;
}
```

```
10x10 lattice

 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1
 1  1  1  1  1  1  1  1  1  1

10x10 lattice

 1 -1 -1  1  1  1 -1  1  1  1
-1  1 -1  1 -1 -1  1 -1 -1  1
 1 -1  1  1 -1 -1  1  1  1  1
 1  1 -1 -1 -1 -1  1 -1 -1  1
 1  1  1 -1 -1 -1 -1 -1  1 -1
-1  1 -1 -1  1  1  1  1  1 -1
-1  1 -1  1 -1 -1 -1 -1 -1 -1
-1  1  1  1  1 -1  1  1  1 -1
 1  1 -1 -1 -1  1 -1 -1 -1  1
-1  1  1 -1  1  1 -1 -1  1  1
```

You can also compute physical values like the energy of the lattice, the magnetization of the lattice, the energy of one particular spin, and the delta of energy of a spin with comparison to the previous state (assuming the spin was flipped and the other were not) :

```cpp
#include "lattice.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
    lattice system(10, 1.);     //  10x10 lattice with temperature 1.

    system.energy();            //  returns the energy
    system.magnetization();     //  returns the magnetization
    system.energy(3, 1);        //  returns the energy of the (3, 1) spin
    system.energy_delta(3, 1);  //  returns the delta of energy of this spin

    return 0;
}
```

Other functions can be found in the [header file](https://github.com/kryzar/Hesiod/blob/master/Program/Program/lattice.hpp) of this class.


### Monte Carlo method and MPI

Each of the following functions is made to work nicely with MPI. That is, it will read *n* simulations at the same time, *n* being the number of cores on your computer. And each expected value is an average of the corresponding value of each distinct simulations running simultaneously with MPI.

The only useful functions are the one which run Monte-Carlo simulations to compute the expected values *<E>*, *<|M|>*, *Cv* and *Khi*. There are two `montecarlo` methods in this class :

#### Fixed temperature

 To use it, just call it with the number of Monte-Carlo cycles you with to create, a folder for the output files and the command line arguments from `main` :

```cpp
#include "lattice.hpp"
#include <string>

int main(int argc, char* argv[])  //  argv must not be const with MPI
{
    lattice system(2, 1.);

    std::string folder = "please/be/indulgent/with/our/results/"; //  folder where the outpt files will appear
    system.montecarlo(5E6, folder, argc, argv); //  5*10^6 Monte-Carlo cycles for this 2x2 lattice

    return 0;
}
```

This first method creates two files :
1. `values(cycle)` contains the expected values as function of the number Monte Carlo cycles, like this :
```
#MC cycles  <E>             <|M|>     variance(E)^2  Cv             Khi                   #accepted configurations
0           -inf            inf       1.28e-05       1.28e-05        3.2e-06              0
250         -8.032          4.016     0.00321264     0.00321264     0.00080316              0
500         -8.016          4.008     0.00641216     0.00641216     0.00160304              0
750       -8.01067        4.00533     0.00961136     0.00961136     0.00240284              0
...
...
...
4999250       -7.98389        3.99464       0.138183       0.138183      0.0184356          10053
4999500       -7.98389        3.99464       0.134997       0.134997       0.017638          10053
4999750       -7.98389        3.99464        0.13181        0.13181      0.0168402          10053
```
2. `energies proba` contains all the energies that appeared in the calculations (each Monte-Carlo cycle consists in swipping one spin, having a new energy for the lattice and see whether we accept or reject this new configuration) with their associated probability of appearance.

Those two files can be used to make plots like the one you can see in this file.

#### Variable temperature and temp-steps

The second method is just a loop of the previous method over an interval of temperatures (the initial temperature is the temperature you initialized the lattice with). It outputs the expected values as functions of the temperature and you simply have to give a final temparature and a temp-step :

```cpp
#include "lattice.hpp"
#include <string>

int main(int argc, char* argv[])  //  argv must not be const with MPI
{
    lattice system(60, 2.);

    std::string folder = "langsam/wozzeck/langsam/"; //  folder where the outpt files will appear
    system.montecarlo(5E6, 3., 0.01, folder, argc, argv); //  5*10^6 cycles per temp-step from 2 to 3 temp-units, with a temp-step of 0.01

    return 0;
}
```

This is what the output file `values(temp)` looks like :

```
2       -7.19972        3.73407        5.80929        1.45232       0.363761         200178
2.01       -9.37298        7.45892       -19.4461       -4.81327       -13.1038         206548
2.02       -16.5375        11.1807       -147.681       -36.1927       -40.1489         208934
2.03        -15.703        14.9031       -114.081       -27.6835       -80.5821         208773
...
...
...
2.97       -301.242        340.922         -88313       -10011.8       -38691.8         592382
2.98       -302.876        344.138       -89003.9       -10022.5       -39297.3         591058
2.99       -300.478        347.343       -87536.9       -9791.49       -39902.9         599835
3       -302.068        350.544       -88726.4       -9858.49       -40510.5         602292
```

The only difference is that the first column is the temperature column.

#### Precautions

- You probably want to run the unit-tests. You simply can add `catch.hpp` (use the version in the repo and not Catch 2 for which this program was not updated to), `unit-tests.hpp` and `unit-tests.cpp` in your project as well, and simply `#include "unit-tests.hpp"`.

- If the output values contain values such as :
```
2       -17796.1    1.4889e+189   1.34078e+148   3.35195e+147           -inf        67206.5
2.005       -35725.5    1.4889e+189   1.34078e+148   3.33526e+147           -inf          74852
2.01       -53477.4    1.4889e+189   1.34078e+148   3.31868e+147           -inf          81042
2.015       -71304.6    1.4889e+189   1.34078e+148   3.30223e+147           -inf        79229.5
2.02       -89050.8    1.4889e+189   1.34078e+148   3.28591e+147           -inf        81443.5
2.025              0              0   1.06281e-314              0    1.2586e+180              0
2.03        -124583    1.4889e+189   1.34078e+148   3.25361e+147           -inf        84168.5
2.035        -142192    1.4889e+189   1.34078e+148   3.23764e+147           
...
```
rerun the simulation. This occures from times to times - we did not manage to solve this - but it more seems to be a randomly appearing bug.
The line with the *0* and the *e-314* is also an outlier. There can be some in the method which goes through many temperatures, but if you choose reasonable temp-step and temp-interval, this should not affect your plots.

- Running the program in release mode will give a significant speed-up.

## License

Be indulgent - this is a school project - and do what you want.
