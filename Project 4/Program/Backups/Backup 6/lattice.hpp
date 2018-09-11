//
//  lattice.hpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#pragma once

#include <fstream>
#include <math.h>
#include "functions.hpp"

#define K 1.3806503E-23 //  boltzmann constant
#define J 3.132E-23     //  J constant

class lattice
{
public:
    
    lattice(void);
    lattice(const unsigned n, const int t);
    lattice(const lattice& other);
    ~lattice();
    
    //  operators
    
    inline int operator()(unsigned row, unsigned col) const;
    inline int operator()(const unsigned* position) const;

    //  modify the latice
    
    void change(const unsigned row, const unsigned col);                   //  change (row, col) to its opposite
    void change(const unsigned* position);                                 //  idem
    void change(const unsigned row, const unsigned col, const int spin);   //  change (row, col) to the input spin
    inline void change_random_spin(void);                                  //  change a random spin
    inline void change_random_spin(unsigned* position);                    //  same, with the position of the random spin
    void positivize(void);                                                 //  sets all the spins to 1
    void negativize(void);                                                 //  sets all the spins to -1
    void randomize(void);                                                  //  randomize all the spins of the lattice
    void montecarlo(const unsigned long mc_cycles, const double final_temp, const double temp_step, const std::string folder);
    
    //  calculations
    
    double energy(void) const;                                                //  energy of the full lattice
    double energy(const unsigned row, const unsigned col) const;              //  energy of (row, col) the spin
    double energy_delta(const unsigned row, const unsigned col) const;        //  delta energie of the spin (assuming it was just flipped)
    double energy_delta(const unsigned* position) const;                      //  idem
    double exp_energy_delta(const unsigned row, const unsigned col) const;    //  returns exp(-beta*delta_energie)
    double exp_energy_delta(const unsigned* position) const;                  //  idem
    double magnetization(void) const;                                         //  magnetization of the lattice
    int neighbors_spins_up(const unsigned row, const unsigned col) const;     //  number of neighbors with up spin
    int neighbors_spins_sum(const unsigned row, const unsigned col) const;    //  sum of the neighbors' spins

    //  outputs
    
    void print(void) const;                 //  prints the lattice to cout
    void print(std::string folder) const;   //  prints to a file in folder
    
private:
    
    //  data
    
    int _dim;           //  the lattice is a _dim x _dim lattice
    double _T;          //  initial temperature
    double* _precalc;   //  array to store the precalculated values of Montecarlo
    int** _config;      //  matrix of all the spins
    
    //  methods
    
    void _metropolis(unsigned* position, double& E, double& M_abs, double& Cv, double& khi, long& accepted_configs);
    void _output(std::ofstream& ofile, const double* averages, const double mc_cycles, const double final_temp, const double temp_step);
    void _output(std::ofstream& out_values, const double* averages, const long mc_cycles, const double Cv, const double khi, const long accepted_configs);
    inline void _update_precalc(void);
    inline void _update_averages(double* averages, const unsigned* position, double& E, double& M_abs);
    
    //  sanity checks

    void _domain_test(const unsigned row, const unsigned col) const;
    void _emptyness_test(void) const;
};

//  operators

inline int lattice::operator()(unsigned row, unsigned col) const
{
    return (_config[row][col]);
}

inline int lattice::operator()(const unsigned* position) const
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

inline void lattice::change_random_spin(unsigned* position)
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

inline void lattice::_update_precalc(void)
{
    double buffer = J / (K * _T);
    
    _precalc[0] = exp(8. * buffer);
    _precalc[1] = exp(4. * buffer);
    _precalc[2] = exp(-4. * buffer);
    _precalc[3] = exp(-8. * buffer);
    _precalc[6] = 1. / (K * _T);
}

inline void lattice::_update_averages(double* averages, const unsigned* position, double& E, double& M_abs)
{
    averages[0] += E;
    averages[1] += E * E;
    averages[2] += abs(M_abs);
    averages[3] += M_abs * M_abs;
}
