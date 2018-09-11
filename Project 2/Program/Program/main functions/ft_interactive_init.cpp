//
//  ft_interactive_init.cpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// libraries
#include <armadillo>
// namespaces
using namespace arma;



mat ft_interactive_init(const long& N, const double& omega, const double& h)
{
    
    mat A(N, N);
    
    for(int i = 0; i < N; i++)
    {
        A(i, i) = 2/(h*h) + omega*omega*h*h*(i*i + 2*i + 1) + 1./((i + 1.)*h);
        if(i != N-1)
        {
            A(i+1, i) = A(i, i+1) = - 1./(h*h);
        }
    }
    
    return (A);
}
