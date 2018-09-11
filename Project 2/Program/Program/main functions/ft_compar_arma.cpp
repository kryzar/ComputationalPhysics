//
//  ft_arma_compar.cpp
//  Program
//
//  Created by Antoine Hugounet on 01/10/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// libraries
#include <fstream>
#include <armadillo>
#include <string>
// headers
#include "time.h"
// namespaces
using namespace arma;
using namespace std;


void ft_arma_compar(ofstream& results, const mat& A)
{
    
    clock_t start;
    clock_t finish;
    double time;
    
    vec eigenvalues;
    
    mat B = A;
    
    start = clock();
    
    eigenvalues = eig_sym(B);
    
    finish = clock();
    time = (double) (finish-start) / ((double) CLOCKS_PER_SEC);
    
    results << "Armadillo time using eig_sym (s): " << time << endl << endl;
}

