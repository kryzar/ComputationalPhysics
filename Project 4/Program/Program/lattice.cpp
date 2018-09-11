//
//  lattice.cpp
//  Program
//
//  Created by Antoine Hugounet on 04/11/2017.
//

#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>        //  std::sort
#include <math.h>
#include "lattice.hpp"
#include "lib/lib.h"
#include "functions.hpp"    //  rng


lattice::lattice(void)
{
    _accepted_configs = 0;
    _Cv = 0.;
    _dim = 0;
    _E = 0.;
    _E_variance = 0.;
    _khi = 0.;
    _mc_cycles = 1; //  avoid "inf" results
    _M = 0.;
    _T = 0.;

    _precalc = nullptr;
    _config = nullptr;
    _averages = nullptr;
}

/////////

lattice::lattice(const unsigned dim, const double T)
{
    _accepted_configs = 0;
    _Cv = 0.;
    _dim = dim;
    _E = 0.;
    _E_variance = 0.;
    _khi = 0.;
    _mc_cycles = 1;
    _M = 0.;
    _T = T;
    double beta = 1. / (K * _T);    //  K is set to 1 but this allows to change it easilly
    
    //  spins initialization
    
    _config = new int*[dim];
    for(int i = 0; i < dim; i++)
    {
        _config[i] = new int[dim];
    }
    
    for(unsigned row = 0; row < _dim; row++)
    {
        for(unsigned col = 0; col < _dim; col++)
        {
            _config[row][col] = 1;
        }
    }
    
    //  delta-energies initializations
    _precalc = new double[7];   //  not very elegant but does the job
    
    _precalc[0] = exp(8. * J * beta);
    _precalc[1] = exp(4. * J * beta);
    _precalc[2] = exp(- 4. * J * beta);
    _precalc[3] = exp(- 8. * J * beta);
    _precalc[4] = 8. * J;
    _precalc[5] = 4. * J;
    _precalc[6] = beta;
    
    _averages = new double[4];
}

/////////

lattice::lattice(const lattice& other)
{
    _accepted_configs = other._accepted_configs;
    _Cv = other._Cv;
    _dim = other._dim;
    _E = other._E;
    _E_variance = other._E_variance;
    _khi = other._khi;
    _mc_cycles = other._mc_cycles;
    _M = other._M;
    _T = other._T;

    
    _config = new int*[_dim];
    for(int i = 0; i < _dim; i++)
    {
        _config[i] = new int[_dim];
    }
    
    for(unsigned row = 0; row < _dim; row++)
    {
        for(unsigned col = 0; col < _dim; col++)
        {
            _config[row][col] = other._config[row][col];
        }
    }
    
    _precalc = new double[7];
    for(int i = 0; i < 7; i++)  _precalc[i] = other._precalc[i];
    
    _averages = new double[4];
    for(int i = 0; i < 4; i++)  _averages[i] = other._averages[i];
    
}

/////////

lattice::~lattice()
{
    for(int i = 0; i < _dim; i++)
    {
        delete[] _config[i];
    }
    
    delete[] _config;
    delete[] _averages;
    delete[] _precalc;
}

//  modify the lattice

void lattice::randomize(void)
{
    _emptyness_test();
    
    for(unsigned row = 0; row < _dim; row++)
    {
        for(unsigned col = 0; col < _dim; col++)
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
    
    for(unsigned row = 0; row < _dim; row++)
    {
        for(unsigned col = 0; col < _dim; col++)
        {
            _config[row][col] = -1;
        }
    }
}

/////////

void lattice::positivize(void)
{
    _emptyness_test();
    
    for(unsigned row = 0; row < _dim; row++)
    {
        for(unsigned col = 0; col < _dim; col++)
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
    _emptyness_test();
    _domain_test(position[0], position[1]);
    
    _config[position[0]][position[1]] *= -1;
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

/////////

void lattice::montecarlo(const unsigned long mc_cycles, const std::string folder, int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    unsigned* position = new unsigned[2];   //  position of the randomly flipped spin
    int rank;
    
    _mc_cycles = mc_cycles;
    std::ofstream out_values(folder + "values(cycle)");
    std::ofstream out_energy(folder + "energies proba");
    std::vector<double> energy_list;    //  vector which contains all the energies appearing, counting multiplicities
    
    if(_T < 1.5)  positivize(); //  ground state spin orientation
    _E = energy();              //  initialize the expected values to be incremented with their variations at each cycle
    _M = magnetization();
    
    for(unsigned long cycle = 0; cycle <= _mc_cycles; cycle++)
    {
        _metropolis(position);      //  makes a move, accept and reject it and computes E and M
        energy_list.push_back(_E);  //  add the energy to the list, again with multiplicities
        _update_averages(position); //  here we update the values of <E>, <|M|>, etc according to what _metropolis did
        _output(out_values, cycle);  //  we use MPI_Reduce here and output one set of averages every 250 cycles
        //  the averages are thus functions of the number of completed cycles
    }

    _energy_proba(out_energy, energy_list, rank);   //  takes the list of energies, sorts it, counts and delete multiplicities
    
    out_values.close();
    out_energy.close();
    delete[] position;
    MPI_Finalize();
}

/////////

void lattice::montecarlo(const unsigned long mc_cycles, const double final_temp, const double temp_step, const std::string folder, int argc, char* argv[])
{
    if(final_temp <= _T || temp_step > (final_temp - _T))
    {
        std::cout << "The final temp must be greater than the initial temp." << std::endl;
        std::cout << "The temp-step must be reasonable." << std::endl;
        exit(6);
    }
    
    //  this functions is just a loop of the above function over temp-steps
    //  see the comments above
    
    MPI_Init(&argc, &argv);
    std::ofstream out_values(folder + "values(temp)");
    unsigned* position = new unsigned[2];   //  position of the randomly flipped spin
    
    _mc_cycles = mc_cycles;
    if(_T < 1.5)    positivize();
    
    while(_T <= final_temp) //  we repeat monte-carlo until the max temp is reached
    {
        _E = energy();
        _M = magnetization();
        _accepted_configs = 0;
        
        for(unsigned long cycle = 0; cycle < _mc_cycles; cycle++)   //  Monte-Carlo
        {
            _metropolis(position);
            _update_averages(position);
        }
        
        //  this time we only output the final average
        // we are interested in the evolution of the quantities with the temperature
        _output(out_values);
        _T += temp_step;
        _update_precalc();
    }
    
    out_values.close();
    delete[] position;
    MPI_Finalize();
}

//  calculations

double lattice::energy(void) const
{
    double energy_lattice = 0.;
    
    for(unsigned row = 0; row < _dim; row++)
    {
        for(unsigned col = 0; col < _dim; col++)
        {
            //  we must not count twice the same spin
            energy_lattice -= (double) _config[row][col] * _config[_periodic(row, _dim, -1)][col] + _config[row][_periodic(col, _dim, -1)];
        }
    }
    
    return (J * energy_lattice);
}

////////

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
    
    //  the conditions represent the boundary conditions
    if(neighbors_sum == 4)          return (- (double) _config[row][col] * _precalc[4]);
    else if(neighbors_sum == 2)     return (- (double) _config[row][col] * _precalc[5]);
    else if(neighbors_sum == 0)     return (0.);
    else if(neighbors_sum == -2)    return ((double) _config[row][col] * _precalc[5]);
    else if(neighbors_sum == -4)    return ((double) _config[row][col] * _precalc[4]);
    
    else
    {
        std::cout << "Error in \"energy_delta\"" << std::endl;
        exit(5);
    }
}


////////

double lattice::energy_delta(const unsigned* position) const
{
    return(energy_delta(position[0], position[1]));
}

////////

double lattice::exp_energy_delta(const unsigned row, const unsigned col) const
{
    _domain_test(row, col);
    _emptyness_test();
    
    int neighbors = neighbors_spins_up(row, col);
    bool positive = _config[row][col] == 1;
    bool negative = _config[row][col] == -1;
    
    //  we use a th precalculations here, they are initialized whe the object is created
    //  those are test to determine in which case we are, depending on the sum of the neighbors
    
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
        return (1.);            //  exp(0) = 1
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
    return (exp_energy_delta(position[0], position[1]));
}

////////

double lattice::magnetization(void) const
{
    double magnet = 0.;
    
    for(unsigned row = 0; row < _dim; row++)
    {
        for(unsigned col = 0; col < _dim; col++)
        {
            magnet += (double) _config[row][col];
        }
    }
    
    return (magnet);
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
    
    if((long) col - 1 < 0)  s1 = _config[row][_dim - 1];
    else                    s1 = _config[row][col-1];
    if((long) row - 1 < 0)  s2 = _config[_dim - 1][col];
    else                    s2 = _config[row - 1][col];
    if(col + 1 > _dim - 1)  s3 = _config[row][0];
    else                    s3 = _config[row][col + 1];
    if(row + 1 > _dim - 1)  s4 = _config[0][col];
    else                    s4 = _config[row + 1][col];
    
    return (s1 + s2 + s3 + s4);
}

//  outputs

void lattice::_energy_proba(std::ofstream& out_energy, std::vector<double>& energy_list, int rank) const
{
    std::vector<std::vector<double>> energy_proba;
    double count;
    unsigned long i;
    unsigned long j;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::sort(energy_list.begin(), energy_list.end());  //  sort the energies
    
    //  the idea is to go through the sorted list of the energies
    //  while the energy is the same as the previous energy, we count the occurences
    //  when the condition is not verifies anymore, we add the energy to the new list
    //  we add its probability as well
    //  we continue this by starting again at the last energy counted
    
    if(rank == 0)
    {
        i = 1;
        
        while(i < energy_list.size())    //  go through the list of energies
        {
            j = i;
            count = 1;
            
            while(abs(energy_list[j] - energy_list[j-1]) == 0.)    // tolerance
            {
                count++;
                j++;
            }
            
            energy_proba.push_back({energy_list[j], count});
            out_energy << energy_proba.back()[0] << setw(15) << energy_proba.back()[1] / _mc_cycles << std::endl;
            i += count;
        }
    }
}

////////

void lattice::_output(std::ofstream& out_values, unsigned long cycle)
{
    if(cycle % 250 == 0)
    {
        int numprocs;
        double* mpi_send = new double[6];       //  <E>, <|M|>, variance(E), Cv, Khi, accepted configs
        double* mpi_receive = new double[6];    //  idem
        
        //  we create a temporary array to send the calculated averages to MPI
        
        mpi_send[0] = _averages[0] / cycle;       //  <E>
        mpi_send[1] = abs(_averages[2]) / cycle;  //  <|M|>
        mpi_send[2] = _E_variance;                //  variance(E)
        mpi_send[3] = _Cv;                        //  Cv
        mpi_send[4] = _khi;                       //  Khi
        mpi_send[5] = _accepted_configs;          //  accepted configs
        
        MPI_Reduce(mpi_send, mpi_receive, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        
        //  we get the sum of the data in mpi_receive
        //  and we make the average with the number of threads
        
        if(mpi_receive[0] != 0)
        {
            out_values << cycle << setw(15);
            out_values << mpi_receive[0] / (double) numprocs  << setw(15);
            out_values << mpi_receive[1] / (double) numprocs  << setw(15);
            out_values << mpi_receive[2] / (double) numprocs  << setw(15);
            out_values << mpi_receive[3] / (double) numprocs  << setw(15);
            out_values << mpi_receive[4] / (double) numprocs  << setw(15);
            out_values << mpi_receive[5] / (double) numprocs  << std::endl;
        }
        
        delete[] mpi_send;
        delete[] mpi_receive;
    }
}

////////

void lattice::_output(std::ofstream& out_values)
{
    //  see above
    
    int numprocs;
    double* mpi_send = new double[6];       //  <E>, <|M|>, variance(E), Cv, Khi, accepted configs
    double* mpi_receive = new double[6];    //  idem
    
    mpi_send[0] = _averages[0] / _mc_cycles;        //  see above
    mpi_send[1] = abs(_averages[2]) / _mc_cycles;
    mpi_send[2] = _E_variance;
    mpi_send[3] = _Cv;
    mpi_send[4] = _khi;
    mpi_send[5] = _accepted_configs;
    
    MPI_Reduce(mpi_send, mpi_receive, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    out_values << _T << setw(15);
    out_values << mpi_receive[0] / (double) numprocs << setw(15);  //  <E>
    out_values << mpi_receive[1] / (double) numprocs << setw(15);  //  <|M|>
    out_values << mpi_receive[2] / (double) numprocs << setw(15);  //  variance(E)
    out_values << mpi_receive[3] / (double) numprocs << setw(15);  //  Cv
    out_values << mpi_receive[4] / (double) numprocs << setw(15);  //  Khi
    out_values << mpi_receive[5] / (double) numprocs << std::endl; //  accepted configs
    
    delete[] mpi_send;
    delete[] mpi_receive;
}

/////////

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
        
        for(unsigned row = 0; row < _dim; row++)
        {
            if(_config[row][0] == 1)
            {
                std::cout << " " ;
            }
            for(unsigned col = 0; col < _dim; col++)
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
        
        for(unsigned row = 0; row < _dim; row++)
        {
            if(_config[row][0] == 1)
            {
                output << " " ;
            }
            for(unsigned col = 0; col < _dim; col++)
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

////////

void lattice::_metropolis(unsigned* position)
{
    change_random_spin(position);   //  flips a random spin and gets its position
    
    //  now we perform the Metropolis tests
    
    if(energy_delta(position) <= 0) //  accept the move
    {
        _E += energy_delta(position);
        _M += 2 * _config[position[0]][position[1]];
        _accepted_configs ++;
    }
    else
    {   //  r <= exp(- beta * delta_energy)
        if(small_random_number() <= exp_energy_delta(position)) //  accept the move
        {
            _E += energy_delta(position);
            _M += 2 * _config[position[0]][position[1]];
            _accepted_configs ++;
        }
        else //  reject the move
        {
            change(position); //  back to the previous state
        }
    }
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
