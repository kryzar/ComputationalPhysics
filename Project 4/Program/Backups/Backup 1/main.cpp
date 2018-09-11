//
//  main.cpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#include <iostream>
#include "lattice.hpp"
#include "functions.hpp"


int main(int argc, const char* argv[])
{
    
    lattice system(20, 1., 1.);
    
    system.print();
    system.randomize();
    system.print();
    
    
    return 0;
}





