#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double** solve_lap_sor(int N, double **bound, double **rho, double w, double tau){
    /* 
    Function for sove laplace equations with SOR method
    N : int, size of grid
    bound : matrix, boundary values
    rho : matrix, source of field
    w : float, overrelaxation parameter
    tau : float, required convergence 
    */
    
    // usefull quantities
    double conv, mean_U0, mean_U, force;
    int iter;
    
    // matrix of the potential
    double **phi;
    phi = (double **)malloc((N+1)*sizeof(double));
    for(int i=0; i<N+1; i++){
        phi[i] = (double*)calloc((N+1), sizeof(double) );
    }
      
    //boundary condition
    for(int i=0; i<=N; i++){
        // attention: the corners could be overwritten
        phi[i][0] = bound[0][i]; // ovest
        phi[i][N] = bound[1][i]; // est
        phi[0][i] = bound[2][i]; // sud
        phi[N][i] = bound[3][i]; // nord
    }

    conv = 10.0;
    mean_U0 = 0;
    iter = 0;
     
    for(int i=0; i<N+1; i++){
        for(int j=0; j<N+1; j++){
            mean_U0 = mean_U0 + phi[i][j];
        }
    }
      
    mean_U0 = mean_U0/((N+1)*(N+1)*1.0);
    
    while(conv > tau){
        for(int i=1; i<N; i++){
            for(int j=1; j<N; j++){
                force = phi[i][j+1] + phi[i][j-1] + phi[i+1][j] + phi[i-1][j];
                force = force + rho[i][j];
                phi[i][j] = w*0.25*force + (1 - w)*phi[i][j];
            }
        }

        mean_U = 0;
        for(int i=0; i<N+1; i++){
            for(int j=0; j<N+1; j++){
                mean_U = mean_U + phi[i][j];
            }
        }

        mean_U = mean_U/((N+1)*(N+1)*1.0);
        conv = fabs(mean_U - mean_U0);
        mean_U0 = mean_U;
        iter = iter + 1;

    }
    
    return phi;
}

int main(void){
    int N = 300;
    double w, tau;
    double **phi, **bound, **rho;

    // potential
    phi = (double **)malloc((N+1)*sizeof(double));
    for(int i=0; i<N+1; i++){
        phi[i] = (double*)calloc((N+1), sizeof(double) );
    }
    
    // boundary condition
    bound = (double **)malloc((4)*sizeof(double));
    for(int i=0; i<4; i++){
        bound[i] = (double*)calloc((N+1), sizeof(double) );
    }
    for(int i=0; i<=N; i++){
        // attention: the corners could be overwritten
        bound[0][i] = 0; // ovest
        bound[1][i] = 0; // est
        bound[2][i] = -2; // sud
        bound[3][i] = 2; // nord
        /*
        if(i <N/4){
            bound[3][i] = 1; // sud
        }
        else if(i > 3*N/4){
            bound[3][i] = 1;
        }
        else{
            bound[3][i] = 0;
        }*/
    }
    
    // Source
    rho = (double **)malloc((N)*sizeof(double));
    for(int i=0; i<N; i++){
        rho[i] = (double*)calloc((N), sizeof(double) );
    }
    for(int i=-1; i<=1; i++){
        for(int j=-1; j<=1; j++){
            rho[7*N/11+i][N/2+j] = -1;
            rho[4*N/11+i][N/2+j] = 1;
        }
    }
       
    w = 1.99;
    tau = 1e-8;

    phi = solve_lap_sor(N, bound, rho, w, tau);

    FILE *fd;
    fd = fopen("lap_c.txt", "w");
    if( fd==NULL ) {
        perror("Erron in opening");
    }

    fprintf(fd, "%d \n", N);

    for(int i=0; i<N+1; i++){
        for(int j=0; j<N+1; j++){
            fprintf(fd, "%.20f \n", phi[i][j]);
        }
    }

    fclose(fd);

    return 0;
}
