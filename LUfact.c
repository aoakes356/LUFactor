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
            printf("%.5f ",A[order[i]][j]);
        }
        printf(" | ");
        printf("\n");
    }
}

void printAugMatrix(double** A, int N, short* order, double* x){
    for(int i = 0; i < N; i++){
        printf("| ");
        for(int j = 0; j < N; j++){
            printf("%.5f ",A[order[i]][j]);
        }
        printf("%f, ",x[i]);
        printf(" | ");
        printf("\n");
    }
    printf("\n");
}

LUfact *LUfactor(int N, const double **A) {
    LUfact *LU = (LUfact*) malloc(sizeof(LUfact));
    LU->N = N;
    LU->LU = createMatrix(N);
    LU->mutate = (short *) malloc(N*sizeof(short));
    double **B = LU->LU;
    float current, k, max = 0;
    int maxIndex;
    short save;
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
    for(int i = 0; i < N-1; i++){ // Column index.
        // Find the maximum pivot in the current column.
        // Swap mutate i and mutate of max index so the pivot row is in the
        // ith index.
        max = 0;
        for(int m = 0; m < N; m++){
            current = B[LU->mutate[m]][i];
            if(abs(current) > max){
                max = current;
                maxIndex = m;
            }
        }
        // Swap here.
        save = LU->mutate[i];
        LU->mutate[i] = LU->mutate[maxIndex];
        LU->mutate[maxIndex] = save;
        // Find the k value here
        for(int j = i+1; j < N; j++){ // Row index.
            k = B[LU->mutate[j]][i]/B[LU->mutate[i]][i];
            B[LU->mutate[j]][i] = k;
            for(int r = i+1; r < N; r++){
                B[LU->mutate[j]][r] -= k*B[LU->mutate[i]][r];
            }
        }
    }
    short a[] = {0,1,2,3,4};
    printf("\n\n"); 
    printMatrix(B,N,LU->mutate); 
    printf("\n\n"); 
    return LU;
}

void LUdestroy(LUfact *fact) {
    destroyMatrix(fact->N,fact->LU);
    free(fact->mutate);
    free(fact);
}

void LUsolve(LUfact *fact, const double *b, double *x) {
    //Solve for LY = B (Lower)
    int N = fact->N;
    double y[N];
    printf("B: ");
    for(int n = 0; n < N; n++){
        printf("%f, ",b[n]);
    }
    printf("\n");
    // Copy b into y.
    for(int n = 0; n < N; n++){
        y[n] = b[fact->mutate[n]];
    }
    for(int j = 0; j < N; j++){
        for(int i = j+1; i < N; i++){
            y[i] -= y[j]*fact->LU[fact->mutate[i]][j];
        }
    }
    printf("Y:  ");
    for(int l = 0; l < N; l++){
        printf("%f,",y[l]);
    }
    printf("\n\n");
    // Copy Y into X
    for(int n = 0; n < N; n++){
        x[n] = y[n];
    }
    for(int i = N-1; i >= 0; i--){
        printf("next = %f/%f \n",x[i],fact->LU[fact->mutate[i]][i]);
        x[i] /= fact->LU[fact->mutate[i]][i];
        for(int j = i-1; j >= 0; j--){
            x[j] -= x[i]*fact->LU[fact->mutate[j]][i];
            printAugMatrix(fact->LU, fact->N, fact->mutate,x);
            fact->LU[fact->mutate[j]][i] = 0;
        }
    }
    printAugMatrix(fact->LU, fact->N, fact->mutate,x);
    printf("X: ");
    for(int n = 0; n < N; n++){
        printf("%f, ",x[n]);
    }
    printf("\n\n");

}
