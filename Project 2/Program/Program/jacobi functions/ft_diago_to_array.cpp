//
//  ft_diago_to_array.cpp
//  Program
//
//  Created by Antoine Hugounet on 29/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// libraries
#include <armadillo>
// namespaces
using namespace arma;
using namespace std;


double* ft_diago_to_array(const mat& A, const long n)
{
 
    double* array = new double[n];
    
    for(long i = 0; i < n; i++)
    {
        array[i] = (double) A(i,i);
    }
    
    sort(array, array + n); // stl function
    
    return (array);
}
