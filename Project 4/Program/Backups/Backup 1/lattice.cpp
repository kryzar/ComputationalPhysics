//
//  lattice.cpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#include <iostream>
#include "lattice.hpp"
#include "lib/lib.h"
#include "functions.hpp"

lattice::lattice(void)
{
    _dim = 0;
    J = 0;
    T = 0;
    _config = nullptr;
}

lattice::lattice(const int n, const int j, const int t)
{
    _dim = n;
    J = j;
    T = t;
    
    _config = new int*[n];
    for(int i = 0; i < n; i++)
    {
        _config[i] = new int[n];
    }
    
    //  spins initialization
    for(int row = 0; row < _dim; row++)
    {
        for(int col = 0; col < _dim; col++)
        {
            _config[row][col] = 1;
        }
    }
}

lattice::lattice(const lattice& other)
{
    _dim = other._dim;
    J = other.J;
    T = other.T;
    
    _config = new int*[_dim];
    for(int i = 0; i < _dim; i++)
    {
        _config[i] = new int[_dim];
    }
    
    for(int row = 0; row < _dim; row++)
    {
        for(int col = 0; col < _dim; col++)
        {
            _config[row][col] = other._config[row][col];
        }
    }
}

lattice::~lattice()
{
    for(int i = 0; i < _dim; i++)
    {
        delete[] _config[i];
    }
    
    delete[] _config;
}

//  methods

int lattice::spin_value(const int row, const int col) const
{
    _emptyness_test();
    
    return (_config[row][col]);
}

void lattice::print(void) const
{
    _emptyness_test();
    
    std::cout << _dim << "x" << _dim << " lattice.\n\n";
    for(int row = 0; row < _dim; row++)
    {
        if(_config[row][0] == 1)
        {
            std::cout << " " ;
        }
        for(int col = 0; col < _dim; col++)
        {
            if((_config[row][col+1] == 1) && (col != _dim -1))
            {
                std::cout << _config[row][col] << "  ";
            }
            else
            {
                std::cout << _config[row][col] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
}

void lattice::randomize(void)
{
    _emptyness_test();
    
    for(int row = 0; row < _dim; row++)
    {
        for(int col = 0; col < _dim; col++)
        {
            if(small_random_number() > 0.5)
            {
                spin_change(row, col);
            }
        }
    }
}

void lattice::spin_change(const int row, const int col)
{
    _emptyness_test();
    
    _config[row][col] = _config[row][col] > 0 ? -1 : 1;
}

void lattice::spin_change(const int row, const int col, const int spin)
{
    _emptyness_test();
    _spin_test(spin);
    
    _config[row][col] = spin;
}

void lattice::_emptyness_test(void) const
{
    if(_dim == 0)
    {
        std::cout << "Empty lattice." << std::endl;
        std::exit(2);
    }
}

void lattice::_spin_test(const int spin)
{
    if(spin != 1 || spin != -1)
    {
        std::cout << "The spin must be 1 or -1" << std::endl;
        std::exit(1);
    }
}
