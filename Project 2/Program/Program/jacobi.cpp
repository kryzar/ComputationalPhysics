//
//  jacobi.cpp
//  Program
//
//  Created by Antoine Hugounet on 18/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//
//  Jacobi's method to compute the eigenvalues of a symmetric matrix


//libraries
#include <fstream>
#include <armadillo>
#include <string>
#include <algorithm>
// headers
#include "jacobi functions/include.hpp"
#include "time.h"
// namespaces
using namespace arma;
using namespace std;


double* jacobi(mat& A, mat& R, const double epsilon, ofstream& stream)
{

    // declarations
    long const n = A.n_cols;
    long k = 0;
    long l = 1;
    long iterations = 0;
    double largest = ft_nondiago_largest(A, n, k, l);
    double c;
    double s;
    double time;
    // the frobenius norm is an invariant of the transformation
    // we use it as a required test
    double frobenius_norm = norm(A, "fro");
    double* eigenvalues = new double[n];
    clock_t start;
    clock_t finish;


    R = eye(n, n);  // reset the eigenmatrix as a security measure
    
    start = clock();
    
    while(fabs(largest) > epsilon && iterations < n*n*n)
    {
        ft_find_trigo(A, c, s, k, l); // finds cos, sin and tan
        ft_rotate(A, n, R, c, s, k, l); // changes the values
        iterations++;
        largest = ft_nondiago_largest(A, n, k, l);
    }
    
    
    finish = clock();
    time = (double) (finish-start) / ((double) CLOCKS_PER_SEC);
    
    // we test the invariant
    
    ft_invariant(A, frobenius_norm, stream);
    
    // sorted dynamic array with the eigenvalues
    
    eigenvalues = ft_diago_to_array(A, n);
    
    // and write the data into a file
    
    ft_cosmetics(stream, A, R, n, time, iterations, eigenvalues);
    
    return (eigenvalues);
}

