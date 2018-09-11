//
//  functions.cpp
//  Program
//
//  Created by Antoine Hugounet on 05/11/2017.
//

#include <chrono>
#include <iostream>
#include <math.h>
#include "lib/lib.h"

double random_number(const int dim)
{
    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    long idum = -now;
    double random = dim * ran3(&idum);
     
    while(random >= dim)
    {
        //  we choose another number so we can use this random as an index
        random = dim * ran3(&idum);
    }
     
    return (random);
}

/////////

double small_random_number(void)
{
    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    long idum = -now;
        
    return (ran3(&idum));
}

/////////

double big_random_number(void)
{
    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    long idum = -now;
    
    return (- std::log(1 - ran3(&idum)));
}
