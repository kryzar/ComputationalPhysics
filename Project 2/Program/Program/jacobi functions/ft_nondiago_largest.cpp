//
//  ft_nondiago_largest.cpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


//libraries
#include <armadillo>
// namespaces
using namespace arma;



double ft_nondiago_largest(mat&A, const long n, long& k, long& l)
{
    
    double largest = 0.;
    
    for(long row = 0; row < n; row++)
    {
        for(long col = row + 1; col < n; col++)
        {
            if((A(row, col) != 0) && (fabs(A(row, col)) > fabs(largest)))
            {
                largest = A(row, col);
                k = row;
                l = col;
            }
        }
    }
    
    return largest;
}


double ft_nondiago_largest(mat& A, const long n) // same without any index
{
    
    double largest = 0.;
    
    for(long row = 0; row < n; row++)
    {
        for(long col = row + 1; col < n; col++)
        {
            if((A(row, col) != 0) && (fabs(A(row, col)) > fabs(largest)))
            {
                largest = A(row, col);
            }
        }
    }
    
    return largest;
}

