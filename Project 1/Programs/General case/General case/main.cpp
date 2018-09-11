//
//  main.cpp
//  General case
//
//  Created by Antoine Hugounet on 08/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// READ ME
// Please feel free to build and run this program!
// The algorithm is in 'gaussian algorithm.cpp' and this files contains a few examples to demonstrate the efficiency, and also some limits of our algorithm.
// On the output file you will find the operation times and the operation time that Armadillo requires for the same operation…


#include "gaussian algorithm.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <string>
#include "time.h"
using namespace arma;


clock_t start, finish, start_LU, finish_LU;
int time_precision {10};
std::string separator="---------------------------------------------------------------------------------";

vec solver(mat& A, vec& w);

// this function is useful not to write a hundred times the same code and make this file cleaner as well
void test(mat& A, vec& w, const int i, std::ofstream& file, const std::string& result="GOOD", const std::string& comment="none.", const bool print="true"){

    mat A_LU; // we have to copy them to compute the LU decomposition using Armadillo and compare the operation times
    vec w_LU;
    vec x_LU;

    if(print){
        file << "A" << i << std::endl << A << std::endl;
        file << "w" << i << std::endl << w << std::endl;
    }

    file << "x" << i << "=solver(A" << i << ", x" << i << ")" << std::endl;
    start=clock();
    file << solver(A, w) << std::endl;
    finish=clock();
    double time = (double) (finish - start)/((double) CLOCKS_PER_SEC);
    file << "Operation time: " << std::setprecision(time_precision) << time << 's' << std::endl;


    start_LU=clock();
    x_LU=solve(A_LU, w_LU);
    finish_LU=clock();
    double time_LU = (double) (finish_LU - start_LU)/((double) CLOCKS_PER_SEC);
    file << "LU decomposition operation time using Armadillo: " << time_LU << 's' << std::endl << std::endl;

    file << "Difference between the two methods (time - time_LU): " << time-time_LU << 's' << std::endl;
    file << "Relative error between the two methods: " << log10(fabs((time_LU-time)/time_LU)) << std::endl << std::endl;


    file << "Result: " << result << std::endl;
    file << "Comments: " << comment << std::endl << std::endl;
    file << separator << std::endl << std::endl;
}



int main(int argc, char* argv[]){

    std::ofstream results; // write the results in a file
    results.open(argv[1]);
    results << "We compute the equation Ai*xi=wi for xi." << std::endl << std::endl << separator << std::endl;

    // here are some tests

    // #1
    mat A1(2,2);
    A1(0,0)=1;
    A1(0,1)=2;
    A1(1,0)=3;
    A1(1,1)=4;

    vec w1(2);
    w1(0)=1;
    w1(1)=7;

    test(A1, w1, 1, results);

    // #3
    mat A2=eye(5,5);

    vec w2(5);
    for(int i=0; i<5; i++) w2(i)=42*i;

    test(A2, w2, 2, results);


    // #4
    mat A3(3,3);
    int value {1};
    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            A3(row, col)=value;
            value++;
        }
    }

    vec w3(3);
    w3(0)=42;
    w3(1)=5;
    w3(2)=-10;

    std::string comment3 {"This matrix is singular and therefore the equation has multiple solutions. The algorithm works for this matrix but it is not reliable.\nIndeed : the w3(0) and w3(1) are the good values, but w3(2) is wrong. You can compute A3*x3 to convince yourself."};

    test(A3, w3, 3, results, "WRONG", comment3, true);


    // #5
    mat A4=eye(10000,10000);

    vec w4(10000);
    for(int i=0; i<10000; i++) w4(i)=(-i+42)*10031;

    std::string comment4 {"However we may have a loss of numerical precision here."};

    test(A4, w4, 4, results, "GOOD", comment4, false);


    results << "Our algorithm can compete Armadillo for tiny n only…" << std::endl;
    results << std::endl << "⊂(￣(ｴ)￣)⊃" << std::endl;



    return 0;
}
