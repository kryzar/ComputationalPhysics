//
//  jacobi.hpp
//  Program
//
//  Created by Antoine Hugounet on 18/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


#pragma once
#include <armadillo>

using namespace arma;
using namespace std;

// prototypes

double* jacobi(mat& A, mat& eigenvec_matrix, const double epsilon, ofstream& stream);
