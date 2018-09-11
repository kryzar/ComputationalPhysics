//
//  include.hpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


#pragma once
#include <armadillo>

using namespace arma;
using namespace std;


//  writes the data to a file including the eigenvalues, times, iterations
void ft_cosmetics(ofstream& stream, const mat& A, const mat& R, const long& n, const double& time, const long& iterations, const double* eigenvalues);

//  exports the diagonal of A to a dynamic sorted array
double* ft_diago_to_array(const mat& A, const long n);

// finds c, s and t given the coordonates of the maximum non-diag element
void ft_find_trigo(mat& A, double& c, double& s, long& k, long& l);

//  checks the stability of the invariant after the Jacobi algorithm computed
//  exits the program if the property is not preserved
void ft_invariant(const mat& A, const double frobenius_norm, ofstream& stream);

//  finds the non-diagonal element of a symmetric matrix
//  overloaded for the unit tests
double ft_nondiago_largest(mat&A, const long n, long& k, long& l);
double ft_nondiago_largest(mat& A, const long n);

//  tests if two given real numbers are equal under any given tolerance
bool ft_numerical_equality(const double x, const double y);

//  once we got our trigonometric values, this function rotates the elements of A
void ft_rotate(mat& A, const long n, mat& R, double& c, double& s, long& k, long& l);
