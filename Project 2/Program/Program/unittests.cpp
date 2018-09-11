//
//  unittests.cpp
//  Program
//
//  Created by Antoine Hugounet on 26/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//


// libraries
#include <armadillo>
#include <string>
#include <fstream>
// headers
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "jacobi.hpp"
#include "jacobi functions/include.hpp"
// namespaces
using namespace std;
using namespace arma;



int run_unittest(int argc, const char* argv[])
{
    
    return Catch::Session().run(argc, argv);
}


TEST_CASE("Unit tests on non_diago_largest()", "[ft_nondiago_largest]"){
    
    mat A = eye(3,3);
    mat B = eye(8,8);
    long k_B;
    long l_B;
    
    A(0,1) = 500420.;
    A(0,2) = -500420.1;
    A(1,2) = -8.;
    
    B(5,7) = 3267887;
    
    REQUIRE(ft_nondiago_largest(A, 3) == -500420.1);
    REQUIRE(ft_nondiago_largest(B, 8) == 3267887);
    
    ft_nondiago_largest(B, 8, k_B, l_B);
    
    REQUIRE(k_B == 5);
    REQUIRE(l_B == 7);
}

TEST_CASE("Unit tests on jacobi()", "[jacobi]")
{
    
    double const epsilon = 1.e-8;
    double* eigenvalues_A = new double[3];
    const string path = "/Users/antoinehugounet/Documents/Scolarité/UiO/FYS3150 - Computational physics/Project 2/Vega/Program/Results/unittests.txt";
    mat A(3, 3);
    mat eigenvectors_A;
    vec eig_A_0;
    vec eig_A_1;
    vec eig_A_2;
    
    ofstream results;
    results.open(path);
    
    A(0,0) = A(1,1) = A(2,2) = 1.;
    A(1,0) = A(0,1) = 2.;
    A(2,0) = A(0,2) = 3.;
    A(2,1) = A(1,2) = 0.;
    
    eigenvalues_A = jacobi(A, eigenvectors_A, epsilon, results);
    
    eig_A_0 = eigenvectors_A.col(0);
    eig_A_1 = eigenvectors_A.col(1);
    eig_A_2 = eigenvectors_A.col(2);
    
    SECTION("Validity of the eigenvalues")
    {
        REQUIRE(ft_numerical_equality(eigenvalues_A[0], 1-sqrt(13)));
        REQUIRE(ft_numerical_equality(eigenvalues_A[1], 1));
        REQUIRE(ft_numerical_equality(eigenvalues_A[2], 1+sqrt(13)));
    }
    
    SECTION("Orthogonality of the eigenvectors")
    {
        REQUIRE(ft_numerical_equality(dot(eig_A_0, eig_A_1), 0));
        REQUIRE(ft_numerical_equality(dot(eig_A_0, eig_A_2), 0));
        REQUIRE(ft_numerical_equality(dot(eig_A_1, eig_A_2), 0));
    }
    
    results.close();
}
