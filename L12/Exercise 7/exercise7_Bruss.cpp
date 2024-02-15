#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

/* our custom libraries */
#include "save_data.h" // Include the custom header file

double sum_vector(vector<double> &x);

void set_matrix_to_zero(int **matrix, int N, int M){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            matrix[i][j] = 0;
        }
    }
}


int main(){
    /*
    ===================================================
    = Exercise 1.8: Gillespie algorithm - Brusselator =
    ===================================================
    */

    // Set up the rnd generator
    int iseed=0; 
    srand(static_cast<unsigned>(iseed));

    // Initial conditions
    int a = 2;
    int b = 5;

    vector<double> Omega = {pow(10, 2), pow(10, 3), pow(10, 4)};

    int t_tot = 1000;
    double dt = pow(10, -3);

    double t = 0.;
    double lambda_C = 0.;
    double u = 0.;
    double tau = 0.; int remain_time;
    int new_state = 0;
    vector <double> p_j = {0., 0., 0., 1.};
    vector<double> rates = {0., 0., 0., 0.};
    int M = 2; int N = t_tot/dt;
    // cout << "N = " << N << endl;
    int **population = (int **)malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        population[i] = (int *)malloc(M * sizeof(int));
    }

    string s = "_brusselator.txt";
    int k = 0;
    vector<int> C = {0, 0};
    
    for (double omega : Omega){
        C[0] = a*int(omega); C[1] =  b*int(omega)/a;
        set_matrix_to_zero(population, N, M);
        k = 0; t = 0;
        population[k][0] = C[0]; population[k][1] = C[1];
        k += 1; t += dt;

        cout << "Omega = " << omega << endl;

        while (t < t_tot){
            rates[0] = a*omega; rates[1] = C[0]; 
            rates[2] = C[0]*(C[0]-1.)*C[1]/(omega*omega); rates[3] = b*C[0];
            lambda_C = sum_vector(rates); // Brusselator equations
            // cout << "lambda_C = " << lambda_C << endl;
            // Extract tau from the exponential distribuation of mean 1/lambda_C
            u = (double)rand() / (RAND_MAX);
            tau = -log(1-u)/lambda_C; 
            // Remain in configuration C for t = tau seconds
            // cout << "tau = " << tau << endl;
            if (tau > dt){
                remain_time = min(int(tau/dt), N-k);
                // cout << "Remain time = " << remain_time << endl;
                for (int i = k; i < k + remain_time; i++){
                    population[i][0] = C[0]; population[i][1] = C[1];
                }
                k += remain_time;
                t += dt*remain_time;
            }
            
            // Then, pick a state j != C
            // For the Brusselator there are 4 reachable states
            p_j[0] = rates[0]/lambda_C; p_j[1] = (rates[0] + rates[1])/lambda_C; 
            p_j[2] = (rates[0] + rates[1] + rates[2])/lambda_C;

            // u = (double)rand() / (RAND_MAX);
            if (u < p_j[0]){
                new_state = 0;
            }
            else if (u < p_j[1] && u > p_j[0]){
                new_state = 1;
            }
            else if (u < p_j[2] && u > p_j[1]){
                new_state = 2;
            }
            else{
                new_state = 3;
            }
            // Change state
            // cout << "State changed in " << new_state << endl;
            if (new_state == 0){
                C[0] += 1;
            }
            else if (new_state == 1){
                C[0] -= 1;
            }
            else if (new_state == 2){
                C[0] += 1; C[1] -= 1;
            }
            else{
                C[0] -= 1; C[1] += 1;
            }   
            // cout << "The new population is C' = " << C[0] << ", " << C[1] << endl; 

            if (k < N){
                population[k][0] = C[0]; population[k][1] = C[1];
            }
            
            t += dt; k += 1;
        }
        // cout << "p_j[2] = " << p_j[2] << endl;
        // cout << "lambda_C = " << lambda_C << endl;
        // cout << "Population C = " << C[0] << " " << C[1] << endl;
        // cout << "Rates :" << rates[0] << " " << rates[1] << " " << rates[2] << " " << rates[3] << " " << endl; 
        saveMatrixToFile_int(population, N, M, "res/" + to_string(int(omega)) + s);

    }

    free(population);

    return 0;
}



double sum_vector(vector<double> &x){
    double res = 0.;
    for(auto val : x){
        res += val;
    }
    return res;
}