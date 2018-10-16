#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUfact.h"

double **createMatrix(int N) {
    double **M = (double **) malloc(N*sizeof(double*));
    for (int i = 0; i < N; i++)
        M[i] = (double*) malloc(N*sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            M[i][j] = (i == j) ? 1.0 : 0.0;
    return M;
}

void destroyMatrix(int N, double **M) {
    for (int i = 0; i < N; i++)
        free(M[i]);
    free(M);
}

LUfact *LUfactor(int N, const double **A) {
    LUfact *LU = (LUfact*) malloc(sizeof(LUfact));
    LU->N = N;
    LU->LU = createMatrix(N);
    LU->mutate = (short *) malloc(N*sizeof(short));
    double **B = createMatrix(N);
    // Clone A into LU
    double **B = LU->LU;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            B[i][j] = B[i][j];

    for (int i = 0; i < N; i++)
        LU->mutate[i] = (short) i;

    // actual factorizing goes here
    short save;
    int b = 0;  // The max possible pivot;
    double k;
    for(int n = 0; n < N; n++){
        // Find the largest pivot in the current collumn.
        for(int k = 0; k < N; k++){
          if(abs(A_[n][k]) > b){
            b = k;
          }  
        }
        // Set the nth value in the mutate array to the index of the largest
        // pivot.
        save = LU->mutate[n];
        LU->mutate[n] = (short)b;
        LU->mutate[b] = save;
        for(int i = n+1; i < N; i++){
            k = 
        }
    }    
    return LU;
}

void LUdestroy(LUfact *fact) {
    destroyMatrix(fact->N,fact->LU);
    free(fact->mutate);
    free(fact);
}

void LUsolve(LUfact *fact, const double *b, double *x) {

}
