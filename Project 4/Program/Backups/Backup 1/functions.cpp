//
//  functions.cpp
//  Program
//
//  Created by Antoine Hugounet on 05/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#include <chrono>
#include <iostream>
#include <math.h>
#include "lib/lib.h"

double small_random_number(void)
{
    //  this is the time since the start of the computer (nano second precision)
    //  it is a long int used as an idum for the ran3 rng
    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    long idum = -now;
        
    return (ran3(&idum));
}

double big_random_number(void)
{
    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    long idum = -now;
    
    return (- std::log(1 - ran3(&idum)));
}
