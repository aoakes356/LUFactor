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

void printMatrix(double** A, int N, short* order){
    for(int i = 0; i < N; i++){
        printf("| ");
        for(int j = 0; j < N; j++){
            printf("%.5f ",A[i][j]);
        }
        printf(" | ");
        printf("\n");
    }
}

LUfact *LUfactor(int N, const double **A) {
    LUfact *LU = (LUfact*) malloc(sizeof(LUfact));
    LU->N = N;
    LU->LU = createMatrix(N);
    LU->mutate = (short *) malloc(N*sizeof(short));
    double **B = createMatrix(N);
    // Clone A into LU
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            B[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < N; i++)
        LU->mutate[i] = (short) i;
    //printf("Initial Matrix \n");
    printMatrix(B,N,LU->mutate);
    // actual factorizing goes here
    short save;
    int b = 0;  // The max possible pivot;
    float max;
    double k;
    for(int n = 0; n < N-1; n++){
        // Find the largest pivot in the current collumn.
        max = 0;
        for(int l = n; l < N; l++){
          if(abs(B[LU->mutate[l]][n]) > max){
            max = B[LU->mutate[l]][n];
            b = l;
          }  
        }
        //printf("Pivot in column %i is %f\n",n,B[b][n]);
        // Set the nth value in the mutate array to the index of the largest
        // pivot.
        //printf("PRE SWAP ---------------\n");
        for(int m = 0; m < N; m++){
            //printf("mutate[%i]: %i\n",m,LU->mutate[m]);
        }
        //printf("Swapping [%i]: %i and [%i]: %i\n",n,LU->mutate[b],b,LU->mutate[n]);
        //printMatrix(B,N,LU->mutate);
        //printf("------------------------\n");
        save = LU->mutate[n];
        LU->mutate[n] = LU->mutate[b];
        LU->mutate[b] = save;
        for(int m = 0; m < N; m++){
            //printf("mutate[%i]: %i\n",m,LU->mutate[m]);
        }
        for(int i = n+1; i < N; i++){
            // Store the k value during each row reduction and put in LU matrix.
            k = B[LU->mutate[i]][n]/B[LU->mutate[n]][n];
            LU->LU[i][n] = k;
            //printf("Current K value: %f \n", k);
            // Reduce
            for(int j = 0; j < N; j++){
                B[LU->mutate[i]][j] -= k*B[LU->mutate[n]][j];
                //printf("Iteration : %i\n",j);
                //printMatrix(B,N,LU->mutate);

            }
            printf("Setting [%i][%i] to %f\n",i,n,k);
            LU->LU[i][n] = k;
            LU->LU[n][i] = B[LU->mutate[n]][i]/B[LU->mutate[n]][n];
            printMatrix(LU->LU, N, LU->mutate);
            printf("----------------------------\n\n");
        }
    }    
    printf("Printing Reduced Matrix B\n");
    printMatrix(B,N, LU->mutate);
    printf("Printing LU matrix \n");
    printMatrix(LU->LU, N, LU->mutate);

    return LU;
}

void LUdestroy(LUfact *fact) {
    destroyMatrix(fact->N,fact->LU);
    free(fact->mutate);
    free(fact);
}

void LUsolve(LUfact *fact, const double *b, double *x) {
    //Solve for Y (Lower)
    double y[fact->N];
    for(int i = 0; i < fact->N-1; i++){
        y[i] = b[i]/fact->LU[i][i];
        fact->LU[i][i] = 1;
        for(int k = i+1; k < fact->N; k++){
            y[k] =b[k] -b[i]*fact->LU[k][i];
            fact->LU[k][i] = 0;
        }
    }
    for(int i = 0; i < fact->N; i++){
        printf("%f, ",y[i]);
    }
    printf("\n");
    //Solve for X (Upper)
    
    for(int i = fact->N-1; i >= 0; i--){
        y[i] /= fact->LU[i][i];
        fact->LU[i][i] = 1;
        for(int k = i-1; k >= 0; k--){
            y[k] -= y[i]*fact->LU[k][i];
            fact->LU[k][i] = 0;
        }
    } 
    for(int i = 0; i < fact->N; i++){
        printf("%f, ",x[i]);
    }
    printf("\n");
}
