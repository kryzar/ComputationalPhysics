//
//  ft_cosmetics.cpp
//  Program
//
//  Created by Antoine Hugounet on 30/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


//libraries
#include <fstream>
#include <string>
// namespaces
using namespace std;


// interactive mode
void ft_cosmetics(ofstream& results, const string mode, const long n, const double rho_max, const double omega)
{

    results << "==== INTERACTIVE CASE ====" << endl << endl;
    results << "Number of meshpoints: " << n+1 << endl;
    results << "rho_max: " << rho_max << endl;
    results << "Omega: " << omega << endl << endl;
}

// non interactive mode
void ft_cosmetics(ofstream& results, const string mode, const long n, const double rho_max)
{
    results << "==== NON INTERACTIVE CASE ====" << endl << endl;
    results << "Number of meshpoints: " << n+1 << endl;
    results << "rho_max: " << rho_max << endl << endl;
}

