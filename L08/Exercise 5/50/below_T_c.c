#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
// #include <cmath>
#include <math.h>

/* our custom libraries */
#include "save_data.h" // Include the custom header file
// #include "RngStream.h"

double compute_energy(int **lattice, double J, int L);
double evaluate_Delta_E(int spin_value, int **lattice, double J, int x, int y, int L);
double sumMatrix(int **matrix, int L);

int main(){
    /*
    ================================
    = Exercise A.6: 2D Ising model =
    ================================
    */

    // Set up the rnd generator
    int iseed=3; //01112023;
    srand(iseed);

    // Ising lattice
    int L = 0; // Linear dimension
    L = 50;
    int **lattice = (int **)malloc(L * sizeof(int *));
    for (int i = 0; i < L; i++) {
        lattice[i] = (int *)malloc(L * sizeof(int));
    }

    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            // To create the random lattice, a random integer is extracted such that σ_i = ±1
            // lattice[i][j] = 2*(rand() % 2) - 1;
            lattice[i][j] = drand48() < 0.5 ? -1:1;
        }
    }



    double J = 1.;

    int x = 0;
    int y = 0;
    int pos = 0;

    int MC_step = L*L;
    int tot_steps = pow(10, 6);
    double Delta_E = 0;
    int spin_value = 0;
    double k_B = 1.;
    double T_c = 2.*J/(k_B*log(1+sqrt(2))); 
    double T = 0.1*T_c;
    double beta = 1./(k_B*T);

    double *Energy_observale = (double*)malloc(tot_steps * sizeof(double));
    // double *Energy_observale_per_step = (double*)malloc(MC_step * sizeof(double));
    // double *Energy_observale_squared = (double*)malloc(tot_steps * sizeof(double));
    // for(int i = 0; i < tot_steps; i++){
    //     // Initialization
    //     Energy_observale_squared[i] = 0;
    // }
    double *Magnetization = (double*)malloc(tot_steps * sizeof(double));
    // double *Magnetization_squared = (double*)malloc(tot_steps * sizeof(double));

    double energy = compute_energy(lattice, J, L);

    for (int t = 0; t < tot_steps; t++){
        for (int mc = 0; mc < MC_step; mc++){
            // Select a random place on the 2D lattice
            pos = rand() % (L*L);   // pos is a random number between 0 and L^2    
            // fprintf(stdout, "Pos is %d", pos);
            y = pos % L;
            x = (pos-y)/L;
            spin_value = lattice[x][y];
            Delta_E = evaluate_Delta_E(spin_value, lattice, J, x, y, L);
            // Metropolis step:
            if(Delta_E < 0 || log(drand48()) < -beta*Delta_E){
                // accept the move
                lattice[x][y] = -spin_value;
                energy += Delta_E;
            }    
        }

        Energy_observale[t] = energy/(L*L);
        Magnetization[t] = sumMatrix(lattice, L)/(L*L);
    }


    
    char s[100] = "res/0_1_T_c_final_lattice_2D_Ising.txt";
    saveMatrixToFile_int(lattice, L, L, s);
    
    char s1[100] = "res/0_1_T_c_magnetization.txt";
    saveResultsToFile(Magnetization, tot_steps, s1);

    char s2[100] = "res/0_1_T_c_energy.txt";
    saveResultsToFile(Energy_observale, tot_steps, s2);

    fprintf(stdout, "Simulation with L = %d and T = %e completed!\n", L, T);

    // free memory
    free(lattice);
    free(Energy_observale);
    free(Magnetization);

    return 0;
}

double compute_energy(int **lattice, double J, int L){
    double E = 0;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            // Interaction with the right neighbour
            E += J*lattice[i][j] * lattice[i][(j+1) % L]; 
            // Interaction with the bottom neighbour
            E += J*lattice[i][j] * lattice[(i+1)%L][j];
        }
    }
    return -E;
}

double evaluate_Delta_E(int spin_value, int **lattice, double J, int x, int y, int L){
    int neighbours_x[4] = {(x-1+L) % L, x, x, (x+1) % L};
    int neighbours_y[4] = {y, (y-1+L) % L, (y+1) % L, y};
    int sum_spins = 0;
    for (int index = 0; index < 4; index++){
        sum_spins += lattice[neighbours_x[index]][neighbours_y[index]];
    }
    return(2*spin_value*J*sum_spins);
}

double sumMatrix(int **matrix, int L) {
    double sum_elements = 0.;
    for(int i = 0; i < L; i++) {
        for(int j = 0; j < L; j++) {
            sum_elements += matrix[i][j];
        }
    }
    return sum_elements;
}
