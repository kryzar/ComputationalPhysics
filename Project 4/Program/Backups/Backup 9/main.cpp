//
//  main.cpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "lattice.hpp"
#include "functions.hpp"
#include "unit-tests.hpp"


int main(int argc, char* argv[])
{
    run_unittest(0, nullptr);
    
    std::string folder = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 4/Hesiod/Results/E/L=40/";
    
    
    lattice system(40, 2.);
    system.montecarlo(1E6, folder, argc, argv);
    
    return 0;
}
