//
//  include.hpp
//  Program
//
//  Created by Antoine Hugounet on 27/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


#pragma once

//  compares the operation time between Armadillo and our algoithm
void ft_arma_compar(ofstream& results, const mat& A);

//  writes the data to a file including the mode, meshpoints, rho_max
//  it is overloaded for each different mode
void ft_cosmetics(ofstream& results, const string mode, const long n, const double rho_max, const double omega);
void ft_cosmetics(ofstream& results, const string mode, const long n, const double rho_max);

//  checks the good number of arguments and that the chosen mode is correct
void ft_precautions(int argc, const char* argv[]);

//  initialize the matrices with the good values for each different mode
mat ft_interactive_init(const long& N, const double& omega, const double& h);
mat ft_noninteractive_init(const long& N, const double& h);
