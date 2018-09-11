//
//  functions.hpp
//  Program
//
//  Created by Antoine Hugounet on 05/11/2017.
//

#pragma once

//  you find here functions that basicly use the ran3 rng from the library
//  they use the high-resolution clock from the C++11 standard library
//  they are just made to make the syntaxe clear

//  creates a random number in [0, dim)
double random_number(const int dim);
//  creates a random number in [0, 1]
double small_random_number(void);
//  creates a random number in [0, infinity)
double big_random_number(void);


