//
//  ft_invariant.cpp
//  Program
//
//  Created by Antoine Hugounet on 29/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// libraries
#include <fstream>
#include <armadillo>
// headers
#include "include.hpp"
// namespaces
using namespace arma;
using namespace std;

void ft_invariant(const mat& A, const double frobenius_norm, ofstream& stream)
{
    
    if(!ft_numerical_equality(frobenius_norm, norm(A, "fro")))
    {
        stream << "ERROR. The Frobenius norm should be invariant." << endl;
        stream << "Results unreliable." << endl;
        exit (2);
    }
}
