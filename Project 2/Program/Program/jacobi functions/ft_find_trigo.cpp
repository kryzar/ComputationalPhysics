//
//  ft_find_trigo.cpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


//libraries
#include <armadillo>
// namespaces
using namespace arma;



void ft_find_trigo(mat& A, double& c, double& s, long& k, long& l)
{
    
    if(A(k,l) != 0)
    {
        double const tau = ((double) A(l,l)-A(k,k))/((double) 2*A(k,l));
        double const t = tau > 0 ? 1./(tau+sqrt(1. + tau*tau)) : -1./(-tau + sqrt(1. + tau*tau));
        c = 1./sqrt(1. + t*t);
        s = c*t;
    }
    else
    {
        c = 1.;
        s = 0.;
    }
}
