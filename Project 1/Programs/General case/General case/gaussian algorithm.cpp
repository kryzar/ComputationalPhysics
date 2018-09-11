//
//  gaussian algorithm.cpp
//  General case
//
//  Created by Antoine Hugounet on 08/09/2017.
//  Copyright © 2017 Hugounet & Villeneuve. All rights reserved.
//

#include "gaussian algorithm.hpp"
#include <armadillo>
using namespace arma;

// we wish to resolve Ax=w, for A a n*n matrix and x and w vectors of size n using the gaussian elimination


vec solver(mat& A, vec& w){
    
    if((w.n_elem!=A.n_cols) || (w.n_elem!=A.n_rows)){
        std::cout << "ERROR. If A is n*n matrix, w must be of size n." << std::endl;
        exit(1);
    }
    
    const long long n=w.n_elem; // size of the vector, long long is an integer which is too big to hold as an int
    
    for(int i=0; i<n; i++){
        if(A(i,i)==0){
            std::cout << "ERROR. A(" << i << "," << i << ")=0" << std::endl;
            exit(2);
        }
    }
    
    // forward substitution, we proceed one row at a time
    
    for(long long row=1; row<n; row++){
        
        double ratio=A(row, row-1)/A(row-1, row-1); // constant ratio on the row
        for(long long i=row; i<n; i++) A(i, row-1)=0; // 0 on the first column of the loop, below the diagonal
        for(long long col=row; col<n; col++) A(row, col)-=ratio*A(row-1, col); // substitutions
        w(row)-=ratio*w(row-1);
        
        if(A(row,row)==0){
            std::cout << "ERROR. Division by 0. " << std::endl;
            exit(3);
        }
    }
    
    // backward substitution
    
    vec x(n);
    
    x(n-1)=w(n-1)/A(n-1,n-1);
    
    for(long long row=n-2; row>=0; --row){
        
        double sum=0.;
        for(long long k=1; k<(n-row); k++) sum+=A(row, n-k)*x(n-k);
        x(row)=(w(row)-sum)/A(row,row);
    }
    
    return x;
}
