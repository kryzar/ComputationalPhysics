//
//  ft_numerical_equality.cpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


//libraries
#include <armadillo>
// namespaces
using namespace arma;



bool ft_numerical_equality(const double x, const double y)
{
    
    double const tolerance = 1.e-4;
    
    return (fabs(x - y) < tolerance);
}
