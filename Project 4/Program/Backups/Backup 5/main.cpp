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
#include "unit-tests.hpp"

using namespace std::chrono;

int main(int argc, const char* argv[])
{
    run_unittest(0, nullptr);
    
    
    lattice test(20, 1.);
    
    test.randomize();
    test.print();
    
    std::cout << test(0, 1) << std::endl;
    
    return 0;
}
