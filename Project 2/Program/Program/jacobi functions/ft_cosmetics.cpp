//
//  ft_cosmetics.cpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


//libraries
#include <fstream>
#include <armadillo>
#include <string>
// namespaces
using namespace arma;
using namespace std;



void ft_cosmetics(ofstream& stream, const mat& A, const mat& R, const long& n, const double& time, const long& iterations, const double* eigenvalues)
{
        
    stream << "==== JACOBI ALGORITHM ====" << endl << endl;
    
    stream << "Matrix size: " << n << 'x' << n << endl;
    stream << "Frobenius norm: OK" << endl;
    stream << "Jacobi operation time (s): " << time << endl;
    stream << "Number of iterations: " << iterations << endl << endl;
    
    stream << "Eigenvalues: " << endl;

    for(int i = 0; i < n; i++)
    {
        stream << eigenvalues[i] << endl;
    }
    
    stream << endl;
    
    if(n < 10)
    {
        stream << "Eigenmatrix (each vector up to be multiplied by a different real constant): " << endl;
        stream << R << endl;
    }
    
    stream << "====/ JACOBI ALGORITHM ====" << endl << endl;
}
