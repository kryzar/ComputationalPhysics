//
//  lattice.cpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include "lattice.hpp"
#include "lib/lib.h"
#include "functions.hpp"

lattice::lattice(void)
{
    _T = 0;
    _dim = 0;
    _precalc = nullptr;
    _config = nullptr;
}

/////////

lattice::lattice(const unsigned n, const int t)
{
    _T = t;
    _dim = n;
    double beta = 1. / (K * _T);
    
    //  spins initialization
    _config = new int*[n];
    for(int i = 0; i < n; i++)
    {
        _config[i] = new int[n];
    }
    
    for(int row = 0; row < _dim; row++)
    {
        for(int col = 0; col < _dim; col++)
        {
            _config[row][col] = 1;
        }
    }
    
    //  delta-energies initializations
    _precalc = new double[6];
    
    _precalc[0] = exp(8. * J * beta);
    _precalc[1] = exp(4. * J * beta);
    _precalc[2] = exp(- 4. * J * beta);
    _precalc[3] = exp(- 8. * J * beta);
    _precalc[4] = 8. * J;
    _precalc[5] = 4 * J;
}

/////////

lattice::lattice(const lattice& other)
{
    _T = other._T;
    _dim = other._dim;
    
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
    
    _precalc = new double[6];
    
    for(int i = 0; i < 6; i++)
    {
        _precalc[i] = other._precalc[i];
    }
    
}

/////////

lattice::~lattice()
{
    for(int i = 0; i < _dim; i++)
    {
        delete[] _config[i];
    }
    
    delete[] _config;
    delete[] _precalc;
}

//  modify the lattice

void lattice::randomize(void)
{
    _emptyness_test();
    
    for(int row = 0; row < _dim; row++)
    {
        for(int col = 0; col < _dim; col++)
        {
            if(small_random_number() > 0.5)
            {
                change(row, col);
            }
        }
    }
}

/////////

void lattice::negativize(void)
{
    _emptyness_test();
    
    for(int row = 0; row < _dim; row++)
    {
        for(int col = 0; col < _dim; col++)
        {
            _config[row][col] = -1;
        }
    }
}

/////////

void lattice::positivize(void)
{
    _emptyness_test();
    
    for(int row = 0; row < _dim; row++)
    {
        for(int col = 0; col < _dim; col++)
        {
            _config[row][col] = 1;
        }
    }
}

/////////

void lattice::change(const unsigned row, const unsigned col)
{
    _emptyness_test();
    _domain_test(row, col);
    
    _config[row][col] *= -1;
}

////////

void lattice::change(const unsigned* position)
{
    int row = position[0];
    int col = position[2];
    
    return (change(row, col));
}

/////////

void lattice::change(const unsigned row, const unsigned col, const int spin)
{
    _emptyness_test();
    _domain_test(row, col);
    
    if(spin != 1 || spin != -1)
    {
        std::cout << "The spin must be 1 or -1" << std::endl;
        std::exit(3);
    }
    
    _config[row][col] = spin;
}

void lattice::montecarlo(const unsigned long mc_cycles, const double final_temp, const double temp_step)
{
    if(final_temp <= _T || temp_step > (final_temp - _T))
    {
        std::cout << "The final temp must be greater than the initial temp." << std::endl;
        std::cout << "The temp step must be reasonable." << std::endl;
        exit(6);
    }
    
    double E;                               // energy
    double M_abs;                           // |magnetic moment| = |M|
    long accepted_configs;
    unsigned* position = new unsigned[2];   // position of the randomly flipped spin
    double* averages = new double[4];       // 0 : <E>, 1 : <E> * <E>, 2 : |M|, 3 : |M| * |M|
    
    for(double temp = _T; temp < final_temp; temp += temp_step)
    {
        E = 0.;
        M_abs = 0.;
        
        //  À INITIALISER, VOIR CODE DE MORTEN
        
        for(long cycle = 0; cycle < mc_cycles; cycle++)
        {
            change_random_spin(position); //  pick a random spin
            
            if(energy_delta(position) <= 0) //  accept the move
            {
                _update_quantities(position, E, M_abs, accepted_configs);
            }
            else
            {   //  r <= exp(-beta*delta_energy)
                if(small_random_number() <= exp_energy_delta(position)) //  accept the move
                {
                    _update_quantities(position, E, M_abs, accepted_configs);
                }
                else //  reject the move
                {
                    change(position); //  back to the previous state
                }
            }
            _update_averages(averages, position, E, M_abs);
        }
        _T += temp_step;
        _update_precalc();
    }
    
    delete[] position;
    delete[] averages;
}

//  calculations

double lattice::energy(const unsigned row, const unsigned col) const
{
    _domain_test(row, col);
    _emptyness_test();
        
    return (- J * (double) _config[row][col] * (double) neighbors_spins_sum(row, col));
}

////////

double lattice::energy_delta(const unsigned row, const unsigned col) const
{
    _domain_test(row, col);
    _emptyness_test();
    
    int neighbors_sum = neighbors_spins_sum(row, col);
    
    if(neighbors_sum == 4)          return (- (double) _config[row][col] * _precalc[4]);
    else if(neighbors_sum == 2)     return (- (double) _config[row][col] * _precalc[5]);
    else if(neighbors_sum == 0)     return (0.);
    else if(neighbors_sum == -2)    return ((double) _config[row][col] * _precalc[5]);
    else if(neighbors_sum == -4)    return ((double) _config[row][col] * _precalc[4]);
    
    else
    {
        std::cout << "Error in \"energy_delta\"" << std::endl;
        exit(5);
        return (-1);
    }
}


////////

double lattice::energy_delta(const unsigned* position) const
{
    int row = position[0];
    int col = position[1];
    
    return(energy_delta(row, col));
}

////////

double lattice::exp_energy_delta(const unsigned row, const unsigned col) const
{
    _domain_test(row, col);
    _emptyness_test();
    
    int neighbors = neighbors_spins_up(row, col);
    bool positive = _config[row][col] == 1;
    bool negative = _config[row][col] == -1;
    
    if((neighbors == 4 && negative) || (neighbors == 0 && positive))
    {
        return (_precalc[3]);   //  exp(-8 * beta * J)
    }
    else if((neighbors == 3 && negative) || (neighbors == 1 && positive))
    {
        return (_precalc[2]);   //  exp(-4 * beta * J)
    }
    else if(neighbors == 2)
    {
        return (1.);            //  exp(0)
    }
    else if((neighbors == 1 && negative) || (neighbors == 3 && positive))
    {
        return (_precalc[1]);   //  exp(4 * beta * J)
    }
    else if((neighbors == 0 && negative) || (neighbors == 4 && positive))
    {
        return (_precalc[0]);   //  exp(8 * beta * J)
    }
    else
    {
        std::cout << "Error in \"exp_energy_delta\"" << std::endl;
        exit(7);
        return (-1);
    }
}

double lattice::exp_energy_delta(const unsigned* position) const
{
    int row = position[0];
    int col = position[1];
    
    return (exp_energy_delta(row, col));
}

////////

int lattice::neighbors_spins_up(const unsigned row, const unsigned col) const
{
    int sum = neighbors_spins_sum(row, col);
    
    if(sum == 4)        return (4);
    else if(sum == 2)   return (3);
    else if(sum == 0)   return (2);
    else if(sum == -2)  return (1);
    else if(sum == -4)  return (0);
    else
    {
        std::cout << "Error in \"neighbors_up_spins\"" << std::endl;
        exit(4);
        return (-1);
    }
}

////////

int lattice::neighbors_spins_sum(const unsigned row, const unsigned col) const
{
    //      s2
    //  s1  X   s3
    //      s4
    //  X = spin(row, col), with periodic conditions
    //  the si are the neighbor's spins
    
    _emptyness_test();
    _domain_test(row, col);
    
    int s1, s2, s3, s4;
    
    if((int) col - 1 < 0)   s1 = _config[row][_dim - 1];
    else                    s1 = _config[row][col-1];
    if((int) row - 1 < 0)   s2 = _config[_dim - 1][col];
    else                    s2 = _config[row - 1][col];
    if(col + 1 > _dim - 1)  s3 = _config[row][0];
    else                    s3 = _config[row][col + 1];
    if(row + 1 > _dim - 1)  s4 = _config[0][col];
    else                    s4 = _config[row + 1][col];
    
    return (s1 + s2 + s3 + s4);
}

//  outputs

void lattice::print(void) const
{
    _emptyness_test();
    
    if(_dim > 40)
    {
        std::cout << "The lattice is too large to be printed." << std::endl;
    }
    else
    {
        std::cout << _dim << "x" << _dim << " lattice\n\n";
        
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
    
}

/////////

void lattice::print(std::string folder) const
{
    _emptyness_test();
    
    std::ofstream output;
    output.open(folder + "lattice");
    
    if(_dim > 70)
    {
        output << "The lattice is too large to be printed." << std::endl;
    }
    else
    {
        output << _dim << "x" << _dim << " lattice\n\n";
        
        for(int row = 0; row < _dim; row++)
        {
            if(_config[row][0] == 1)
            {
                output << " " ;
            }
            for(int col = 0; col < _dim; col++)
            {
                if((_config[row][col+1] == 1) && (col != _dim -1))
                {
                    output << _config[row][col] << "  ";
                }
                else
                {
                    output << _config[row][col] << " ";
                }
            }
            output << std::endl;
        }
        output << std::endl;
    }
    output.close();
}


//  sanity checks

void lattice::_domain_test(const unsigned row, const unsigned col) const
{
    if(row >= _dim || col >= _dim)
    {
        std::cout << "Wrong domain index." << std::endl;
        std::exit(1);
    }
}

/////////

void lattice::_emptyness_test(void) const
{
    if(_dim == 0)
    {
        std::cout << "Empty lattice." << std::endl;
        std::exit(2);
    }
}
