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
    //  this is the time since the start of the computer (nano second precision)
    //  it is a long int used as an idum for the ran3() rng
    static auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    static long idum = -now;
    double random = dim * ran3(&idum);  //  do not declare it as static
     
    while(random == dim)
    {
        //  in the case where ran3(&idum) == 1, very unlikely…
        //  we choose another number so we can use this random as an index
        random = dim * ran3(&idum);
    }
     
    return (random);
}

/////////

double small_random_number(void)
{
    static auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    static long idum = -now;
        
    return (ran3(&idum));
}

/////////

double big_random_number(void)
{
    static auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    static long idum = -now;
    
    return (- std::log(1 - ran3(&idum)));
}
