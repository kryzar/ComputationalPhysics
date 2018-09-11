//
//  lattice.hpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//

#pragma once

#include <fstream>
#include <math.h>
#include "functions.hpp"

#define K 1    //  boltzmann constant, 1.39951E-23
#define J 1.   //  J constant, 1.39951E-23

class lattice
{
public:
    
    lattice(void);
    lattice(const unsigned dim, const double T);
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
    void montecarlo(const unsigned long mc_cycles, const std::string folder, int argc, char* argv[]);
    void montecarlo(const unsigned long mc_cycles, const double final_temp, const double temp_step, const std::string folder, int argc, char* argv[]);
    
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
    
    int _dim;                   //  lattice dimension
    double _Cv;                 //  specific heat
    double _E;                  //  energy
    double _E_variance;         //  variance of E
    double _khi;                //  susceptibility
    double _M;                  //  |magnetic moment| = |M|
    double _T;                  //  initial temperature
    unsigned long _mc_cycles;   //  monte carlo cycles
    unsigned long _accepted_configs;     //  accepted configurations in Monte Carlo
    double* _averages;          //  <E>, <E*E>, <|M|>, <|M|*|M|>
    double* _precalc;           //  array to store the precalculated values of Montecarlo
    int** _config;              //  matrix of all the spins
    
    void _domain_test(const unsigned row, const unsigned col) const;
    void _emptyness_test(void) const;
    void _output(std::ofstream& out_values, double* mpi_send, double* mpi_receive);
    void _output(std::ofstream& out_values, double* mpi_send, double* mpi_receive, unsigned long cycle);
    void _metropolis(unsigned* position);
    inline void _update_precalc(void);
    inline void _update_averages(const unsigned* position);
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
    change(random_number(_dim), random_number(_dim));
}

inline void lattice::change_random_spin(unsigned* position)
{
    position[0] = (unsigned) random_number(_dim);
    position[1] = (unsigned) random_number(_dim);
    
    change(position[0], position[1]);
}

inline void lattice::_update_precalc(void)
{
    double buffer = J / (K * _T);
    
    _precalc[0] = exp(8. * buffer);
    _precalc[1] = exp(4. * buffer);
    _precalc[2] = exp(-4. * buffer);
    _precalc[3] = exp(-8. * buffer);
    _precalc[6] = 1. / (K * _T);    //  beta
}

inline void lattice::_update_averages(const unsigned* position)
{
    _averages[0] += _E;
    _averages[1] += _E * _E;
    _averages[2] += _M;
    _averages[3] += _M * _M;
    
    _E_variance = (_averages[1] / _mc_cycles) - (_averages[0] / _mc_cycles)*(_averages[0] / _mc_cycles);
    _Cv = (_precalc[6] / _T) * _E_variance;
    _khi = _precalc[6] * ((_averages[3] / _mc_cycles) - (_averages[2] / _mc_cycles)*(_averages[2] / _mc_cycles));
}
