//
//  lattice.hpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#pragma once

class lattice
{
public:
    lattice(void);
    lattice(const int n, const int j, const int t);
    lattice(const lattice& other);
    ~lattice();
    
    //  data
    
    double J;
    double T;
    
    //  method
    
    int spin_value(const int row, const int col) const; //  returns the (row, col) spin value
    void print(void) const; //  prints the lattice to cout
    void randomize(void);   //  randomize all the spins of the lattice
    void spin_change(const int row, const int col); //  change the spin (row, col) to its opposite
    void spin_change(const int row, const int col, const int spin); //  change the spin (row, col) to the new input spin
    
private:
    int _dim;
    int** _config;
    void _emptyness_test(void) const;
    void _spin_test(const int spin);
};
