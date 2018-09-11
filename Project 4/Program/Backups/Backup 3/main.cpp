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
    
    return 0;
}
