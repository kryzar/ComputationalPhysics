//
//  lattice.hpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#pragma once
#include <math.h>
#include "functions.hpp"

#define K 1.3806503E-23 //  boltzmann constant
#define J 3.132E-23     //  J constant

class lattice
{
public:
    
    lattice(void);
    lattice(const int n, const int t);
    lattice(const lattice& other);
    ~lattice();
    
    //  operators
    
    inline int operator()(int row, int col) const;
    inline int operator()(int* position) const;

    //  modify the latice
    
    void change(const int row, const int col);                   //  change (row, col) to its opposite
    void change(int* position);                                   //  idem
    void change(const int row, const int col, const int spin);  //  change (row, col) to the input spin
    inline void change_random_spin(void);                        //  change a random spin
    inline void change_random_spin(int* position);               //  same, with the position of the random spin
    void positivize(void);                                        //  sets all the spins to 1
    void negativize(void);                                        //  sets all the spins to -1
    void randomize(void);                                         //  randomize all the spins of the lattice
    void montecarlo(const unsigned long mc_cycles, const double final_temp, const double temp_step);
    
    //  calculations
    
    double energy(const int row, const int col) const;              //  energie of (row, col) the spin
    double energy_delta(const int row, const int col) const;        //  delta energie of the spin (assuming it was just flipped)
    double energy_delta(int* position) const;                        //  idem
    double exp_energy_delta(const int row, const int col) const;    //  returns exp(-beta*delta_energie)
    double exp_energy_delta(int* position) const;                    //  idem
    int neighbors_spins_up(const int row, const int col) const;     //  number of neighbors with up spin
    int neighbors_spins_sum(const int row, const int col) const;    //  sum of the neighbors' spins

    //  outputs
    
    void print(void) const;                                 //  prints the lattice to cout
    void print(std::string folder) const;                   //  prints to a file in folder
    
private:
    
    //  data
    
    int _dim;               //  the lattice is a _dim x _dim lattice
    double _T;              //  initial temperature
    double* _precalc;       //  array to store the precalculated values of Montecarlo
    int** _config;          //  matrix of all the spins
    
    //  methods
    
    inline int _spin(const int row, const int col) const;
    inline int _spin(int* position) const;
    inline void _update_precalc(void);
    inline void _update_quantities(int* position, double& E, double& M_abs, long& accepted_configs);
    inline void _update_averages(double* averages, int* position, double& E, double& M_abs);
    
    //  sanity checks

    void _domain_test(const int row, const int col) const;
    void _emptyness_test(void) const;
};

//  operators

inline int lattice::operator()(int row, int col) const
{
    return (_config[row][col]);
}

inline int lattice::operator()(int* position) const
{
    unsigned row = position[0];
    unsigned col = position[1];
    
    return (_config[row][col]);
}


//  inline methods


inline void lattice::change_random_spin(void)
{
    _emptyness_test();
    
    change((int) random_number(_dim), (int) random_number(_dim));
}

inline void lattice::change_random_spin(int* position)
{
    _emptyness_test();
    
    static int x = (int) random_number(_dim);
    long n = 0;
    for(long i = 0; i < 100000000; i++)
    {
        n++;
    }
    static int y = (int) random_number(_dim);
    position[0] = x;
    position[1] = y;
    change(x, y);
}

inline int lattice::_spin(const int row, const int col) const
{
    _emptyness_test();
    _domain_test(row, col);
    
    return (_config[row][col]);
}

inline int lattice::_spin(int* position) const
{
    int row = position[0];
    int col = position[1];
    
    return (_config[row][col]);
}

inline void lattice::_update_precalc(void)
{
    double buffer = J / (K * _T);
    
    _precalc[0] = exp(8. * buffer);
    _precalc[1] = exp(4. * buffer);
    _precalc[2] = exp(-4. * buffer);
    _precalc[3] = exp(-8. * buffer);
}

inline void lattice::_update_quantities(int* position, double& E, double& M_abs, long& accepted_configs)
{
    E += energy_delta(position);
    M_abs += 2 * _spin(position);
    accepted_configs ++;
}

void lattice::_update_averages(double* averages, int* position, double& E, double& M_abs)
{
    averages[0] += E;
    averages[1] += E * E;
    averages[2] += M_abs;
    averages[3] += M_abs * M_abs;
}
