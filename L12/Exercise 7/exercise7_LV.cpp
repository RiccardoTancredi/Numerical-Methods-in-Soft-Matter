#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

/* our custom libraries */
#include "save_data.h" // Include the custom header file



int main(){
    /*
    =========================================================
    = Exercise 1.8: Gillespie algorithm -- Lotka - Volterra =
    =========================================================
    */

    // Set up the rnd generator
    int iseed=3; 
    srand(static_cast<unsigned>(iseed));

    // Initial conditions
    vector<int> C = {500, 75};

    double k_1 = 3.;
    double k_2 = 1./100;
    double k_3 = 5.;

    int t_tot = 10;
    double dt = pow(10, -4);

    double t = 0.;
    double lambda_C = 0.;
    double u = 0.;
    double tau = 0.; int remain_time;
    int new_state = 0;
    vector <double> p_j = {0., 0., 1.};
    vector<double> rates = {0., 0., 0.};
    int M = 2; int N = t_tot/dt;
    cout << "N = " << N << endl;
    int **population = (int **)malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        population[i] = (int *)malloc(M * sizeof(int));
    }

    string s = "population_no_stat.txt";
    int k = 0;

    population[k][0] = C[0]; population[k][1] = C[1];
    k += 1; t += dt;
    
    while (t < t_tot){
        lambda_C = k_1*C[0] + k_2*C[0]*C[1] + k_3*C[1]; // Lotka - Volterra equations
        rates[0] = k_1*C[0]; rates[1] = k_2*C[0]*C[1]; rates[2] = k_3*C[1];

        // Extract tau from the exponential distribuation of mean 1/lambda_C
        u = (double)rand() / (RAND_MAX);
        tau = -log(1-u)/lambda_C; 
        // Remain in configuration C for t = tau seconds
        cout << "tau = " << tau << endl;
        if (tau > dt){
            remain_time = min(int(tau/dt), N-k);
            cout << "Remain time = " << remain_time << endl;
            for (int i = k; i < k + remain_time; i++){
                population[i][0] = C[0]; population[i][1] = C[1];
            }
            k += remain_time;
            t += dt*remain_time;
        }
        
        // Then, pick a state j != C
        // For the Lotka - Volterra there are 3 reachable states
        p_j[0] = rates[0]/lambda_C; p_j[1] = p_j[0] + rates[1]/lambda_C; 

        u = (double)rand() / (RAND_MAX);
        if (u < p_j[0]){
            new_state = 0;
        }
        else if (u < p_j[1] && u > p_j[0]){
            new_state = 1;
        }
        else{
            new_state = 2;
        }
        // Change state
        cout << "State changed in " << new_state << endl;
        if (new_state == 0){
            C[0] += 1;
        }
        else if (new_state == 1){
            C[0] -= 1; C[1] += 1;
        }
        else{
            C[1] -= 1;
        }   
        cout << "The new population is C' = " << C[0] << ", " << C[1] << endl; 

        if (k < N){
            population[k][0] = C[0]; population[k][1] = C[1];
        }
        
        t += dt; k += 1;
    }
    
    saveMatrixToFile_int(population, N, M, s);

    free(population);

    return 0;
}

