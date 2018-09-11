//
//  functions.cpp
//  General case
//
//  Created by Antoine Hugounet on 01/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//

#include "functions.hpp"
#include <cmath>

double f(double x){
    return (100.*exp(-10.*x));
}

double u(double x){
    return (1.-(1.-exp(-10.))*x-exp(-10.*x));
}
