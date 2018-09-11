//
//  ft_rotate.cpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


//libraries
#include <armadillo>
// namespaces
using namespace arma;



void ft_rotate(mat& A, const long n, mat& R, double& c, double& s, long& k, long& l)
{
        
    double A_kk = A(k,k);
    double A_ll = A(l,l);
    double A_rowk;
    double A_rowl;
    double R_rowk;
    double R_rowl;
        
    A(k,k) = (double) A_kk*c*c - 2*A(k,l)*c*s + A_ll*s*s;
    A(l,l) = (double) A_ll*c*c + 2*A(k,l)*c*s + A_kk*s*s;
    A(k,l) = 0.;
    A(l,k) = 0.;
    
    for(long row = 0; row < n; row++)
    {
        if(row != k && row != l)
        {
            A_rowk = A(row,k);
            A_rowl = A(row,l);
            A(row,k) = A(k,row) = (double) A_rowk*c - A_rowl*s;
            A(row,l) = A(l,row) = (double) A_rowl*c + A_rowk*s;
            
        }
        // set up the eigenvector
        R_rowk = R(row, k);
        R_rowl = R(row, l);
        R(row, k) = (double) R_rowk*c - R_rowl*s;
        R(row, l) = (double) R_rowl*c + R_rowk*s;;
    }
}
