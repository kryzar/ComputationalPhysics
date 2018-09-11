//
//  tests.cpp
//  Program
//
//  Created by Antoine Hugounet on 26/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//
//  Some basic tests to be computed in main() if you want so


// libraries
#include <armadillo>
#include <fstream>
// headers
#include "jacobi.hpp"
// namespaces
using namespace std;
using namespace arma;


void ft_test_I(const double epsilon, ofstream& stream)
{
    
    mat I = eye(52, 52);
    mat eigen_I;
    
    jacobi(I, eigen_I, epsilon, stream);
}

void ft_test_A(const double epsilon, ofstream& stream)
{
    
    mat A(2,2);
    mat eigen_A;
    A(0, 0) = 3;
    A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = 3;
    
    jacobi(A, eigen_A, epsilon, stream);
}

void ft_test_B(const double epsilon, ofstream& stream)
{
    
    mat B(2,2);
    mat eigen_B;
    B(0, 0) = 1;
    B(0, 1) = 2;
    B(1, 0) = 2;
    B(1, 1) = 1;
    
    jacobi(B, eigen_B, epsilon, stream);
}

void ft_test_C(const double epsilon, ofstream& stream)
{
    
    mat C(3,3);
    mat eigen_C;
    C(0, 0) = C(1, 1) = C(2, 2) = 1.;
    C(1, 0) = C(0, 1) = 2.;
    C(2, 0) = C(0, 2) = 3.;
    C(2, 1) = C(1, 2) = 0.;
    
    jacobi(C, eigen_C, epsilon, stream);
}

void ft_tests(void)
{
    
    const string path = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Program/tests";
    ofstream stream;
    stream.open(path);
    
    const double epsilon = 1.e-8;
    
    ft_test_I(epsilon, stream);
    ft_test_A(epsilon, stream);
    ft_test_B(epsilon, stream);
    ft_test_C(epsilon, stream);

    stream.close();
}
