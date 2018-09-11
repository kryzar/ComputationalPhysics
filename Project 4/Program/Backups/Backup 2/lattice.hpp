//
//  lattice.hpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#pragma once
#include "functions.hpp"

class lattice
{
public:
    
    lattice(void);
    lattice(const int n, const int j, const int t);
    lattice(const lattice& other);
    ~lattice();
    
    //  getters
    
    int accepted_config(void) const;                        //  number of accepted configurations in MC
    inline int value(const int row, const int col) const;   //  returns the (row, col) spin

    
    //  modify the latice
    
    void change(const int row, const int col);                  //  change (row, col) to its opposite
    void change(const int row, const int col, const int spin);  //  change (row, col) to the input spin
    inline void change_random_spin(void);                       //  change a random spin
    inline void change_random_spin(int* position);              //  same, with the position of the random spin
    void positivize(void);                                      //  sets all the spins to 1
    void negativize(void);                                      //  sets all the spins to -1
    void randomize(void);                                       //  randomize all the spins of the lattice
    
    //  calculations
    
    double energie(const int row, const int col) const;             //  energie of (row, col) the spin
    double energie_delta(const int row, const int col) const;       //  delta energie of the spin (assuming it was just flipped)
    int neighbors_spins_up(const int row, const int col) const;     //  number of neighbors with up spin
    int neighbors_spins_sum(const int row, const int col) const;    //  sum of the neighbors' spins

    //  outputs
    
    void print(void) const;                                 //  prints the lattice to cout
    void print(std::string folder) const;                   //  prints to a file in folder
    
private:
    
    //  data
    
    int _dim;               //  the lattice is a _dim x _dim lattice
    int _accepted_config;   //  number of accepted configurations in a MC cycle
    double _J;
    double _T;
    double* _precalc;       //  array to store the precalculated values
    int** _config;          //  matrix of all the spins
    
    //  sanity checks

    void _domain_test(const int row, const int col) const;
    void _emptyness_test(void) const;
};


//  inline methods

inline void lattice::change_random_spin(void)
{
    _emptyness_test();
    
    change((int) random_number(_dim), (int) random_number(_dim));
}

inline void lattice::change_random_spin(int* position)
{
    _emptyness_test();
    
    int x = (int) random_number(_dim);
    int y = (int) random_number(_dim);
    position[0] = x;
    position[1] = y;
    change(x, y);
}

inline int lattice::value(const int row, const int col) const
{
    _emptyness_test();
    _domain_test(row, col);
    
    return (_config[row][col]);
}

