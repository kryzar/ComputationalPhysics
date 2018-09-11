//
//  ft_precautions.cpp
//  Program
//
//  Created by Antoine Hugounet on 29/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// headers
#include <iostream>
#include "unittests.hpp"
// namespaces
using namespace std;



void ft_precautions(int argc, const char* argv[])
{
    
    string const mode = argv[1]; // interactive or non-interactive mode
    
    if(run_unittest(0, nullptr) != 0) // avoid collision with main arguments
    {
        cout << "Unittests failed." << endl;
        exit(1);
    }
    
    if(mode != "interactive" && mode != "noninteractive")
    {
        cout << "Choose the interactive or noninteractive mode." << endl;
        exit(10);
        
    }
    else if(mode == "interactive" && argc != 5)
    {
        cout << "5 cmd. line arguments required in interactive mode." << endl;
        cout << "name, mode, number of mesh points, rho_max, omega" << endl;
        exit(11);
    }
    else if(mode == "noninteractive" && argc != 4)
    {
        cout << "4 cmd. line arguments required in interactive mode." << endl;
        cout << "name, mode, number of mesh points, rho_max" << endl;
        exit(12);
    }
    

}
