//
//  main.cpp
//  General case
//
//  Created by Antoine Hugounet and Ethel Villeneuve on 01/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//

#include <iostream>
#include <cmath>        // math functions
#include <fstream>      // write in files
#include <iomanip>      // write in files
#include "functions.hpp" // contains function f and the exact u, for more readability
#include <string>
#include "time.h"


double f(double x);
double u(double x);


int main(int argc, char* argv[]){

    // FIRST PART, THE GAUSSIAN ELIMINATION TO OBTAIN OUR VECTOR V

    if(argc!=7){
        std::cout << "ERROR. This program must work with 7 command-line arguments for main." << std::endl;
        exit(1);
    }

    int n = atoi(argv[1]);

    double h = 1/(double(n)+1);
    double hSquare = 1/((double)n*n+2*n+1); // so the computer does not need to calculate this value a billion times

    double* g = new double[n]; // dynamic array for b tilde
    for(long i = 0; i<n; i++) g[i]=hSquare*f((double(i)+1)*h);

    // WARNING!
    // g[0] in the code does not correspond to g_0 in the reality.
    // it is the first element of the array g, but it contains the value of b_1 tilde
    // g[i] in the program in fact corresponds to g_(i+1) in the reality

    double* a = new double[n-1]; // lower diagonal
    double* b = new double[n]; // main diagonal
    double* c = new double[n-1]; // upper diagonal
    for(long i=0; i < n-1; i++) a[i] = atof(argv[2]);
    for(long i=0; i < n; i++) b[i] = atof(argv[3]);
    for(long i=0; i < n-1; i++) c[i] = atof(argv[4]);

    clock_t start, finish;

    // this step consists in making A an upper-triangular matrix by forward substitution…

    start=clock();
    for(long i=1; i<n; i++){g
        b[i]-=(a[i-1]/b[i-1])*c[i-1];
        g[i]-=(a[i-1]/b[i-1])*g[i-1];
        // no need to change the values of a, even if they change mathematically
    }

    // …then backward substitution

    double* v = new double[n]; // creation of the dynamic array containing the v_i, idem, v[i] in the program is v_(i+1) in the reality

    v[n-1]=g[n-1]/b[n-1];
    for(long i=n-2; i>=0; --i) v[i]=(g[i]-c[i]*v[i+1])/b[i];

    // free the memory
    delete[] a;
    delete[] b;
    delete[] c;
    finish=clock();

    // SECOND PART, EVALUATION OF THE APPROXIMATION

    std::ofstream results; // open a stream for the results file
    results.open(argv[5]);
    std::ofstream data; // open a stream for the data file, a quick preview of the results file which is heavy for large n
    data.open(argv[6]);

    std::string separator="-------------------------------------------------------------------------------";

    results << "Gaussian elimination results for a tridiagonal matrix of size n=" << n << ".\nWe have x(0)=0 and x_(n+1)=1." << std::endl << std::endl;
    results << separator << std::endl;
    results << "x value                     num. solution       exact solution      relative error" << std::endl;
    results << separator << std::endl;

    int width=20;
    int precision=10;

    double maxerror=log10(fabs((v[0]-u(1./((double)n+1))/v[0])));
    long mesh=0;

    for(long i=0; i<n; i++){

        double xi=(double(i)+1.)/(double(n)+1.);
        double error=log10(fabs(((v[i]-u(xi))/v[i])));
        if(error<maxerror){
            maxerror=error;
            mesh=i;
        }

        results << std::setprecision(precision) << "x" << i+1;
        results << std::setw(width) << std::setprecision(precision) << xi;                        // xi
        results << std::setw(width) << std::setprecision(precision) << v[i];                       // num. solution
        results << std::setw(width) << std::setprecision(precision) << u(xi);                      // exact solution
        results << std::setw(width) << std::setprecision(precision) << error << std::endl;         // relative error
    }

    results << std::endl << "Maximum relative error: " << maxerror << " for i=" << mesh << "." << std::endl;
    results << "Operation time: " << (double) (finish-start)/((double) CLOCKS_PER_SEC) << 's' << std::endl;
    results << std::endl << std::endl << "ヽ(°〇°)ﾉ " << std::endl;

    data << "n=" << n << std::endl << std::endl;
    data << "Maximum relative error: " << maxerror << " for i=" << mesh << "." << std::endl;
    data << "Operation time: " << (double) (finish-start)/((double) CLOCKS_PER_SEC) << 's' << std::endl;
    data << std::endl << std::endl << "(#ಠQಠ#) " << std::endl;

    results.close();
    data.close();

    delete[] g;
    delete[] v;

    return 0;
}
