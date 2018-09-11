//
//  main.cpp
//  Program
//
//  Created by Antoine Hugounet on 18/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// libraries
#include <fstream>
#include <armadillo>
#include <string>
// headers
#include "jacobi.hpp"
#include "ft_tests.hpp" // some basic tests, use ft_tests() here to see them
#include "main functions/include.hpp"
#include "unittests.hpp"
// namespaces
using namespace arma;
using namespace std;


int main(int argc, const char* argv[])
{
    
    // argv[1]: mode [string] "interactive" or "noninteractive"
    // argv[2]: n [int] number of meshpoints
    // argv[3]: rho_max [double]
    // argv[4]: omega [double] only in interactive mode
    
    ft_precautions(argc, argv);
        
    long const n = atol(argv[2]) - 1; // we want n from 1 to N-1
    double const rho_max = atof(argv[3]);
    double const h = rho_max / n;
    double const epsilon = 1.e-8;
    string const mode = argv[1]; // interactive or non-interactive mode
    string const path = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Program/Results/last results.txt";
    mat A(n, n);
    mat R(n, n);

    
    ofstream results;
    results.open(path);
    
    if(mode == "interactive")
    {
        const double omega = atof(argv[4]);

        ft_cosmetics(results, mode, n, rho_max, omega); // a function which writes data to a file
        A = ft_interactive_init(n, omega, h);           // initialisation of A
        //  ft_tests();                                 // run some tests
        //  ft_arma_compar(results, A);                 // compare with Armadillo
        jacobi(A, R, epsilon, results);                 // jacobi algorithm
    }
    else
    {
        ft_cosmetics(results, mode, n, rho_max);
        A = ft_noninteractive_init(n, h);
        jacobi(A, R, epsilon, results);
    }
    
    results.close();

    return (0);
    
}
