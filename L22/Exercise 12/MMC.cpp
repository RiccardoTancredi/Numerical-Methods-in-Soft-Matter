#include <iostream>
#include <set>
#include <unordered_set>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

/* our custom libraries */
#include "save_data.h" // Include the custom header file

double compute_energy(int **lattice, double J, int L);
double evaluate_Delta_E(int spin_value, int **lattice, double J, int x, int y, int L);
double sumMatrix(int **matrix, int L);


int** generateRandomLattice(int L) {
    int** lattice = new int*[L];
    for (int i = 0; i < L; i++) {
        lattice[i] = new int[L];
    }

    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            // To create the random lattice, a random integer is extracted such that σ_i = ±1
            // lattice[i][j] = 2 * (rand() % 2) - 1;
            lattice[i][j] = (double)rand() / (RAND_MAX) < 0.5 ? -1 : 1;
        }
    }

    return lattice;
}


// Function to perform the Metropolis update
template <typename LatticeType>
void metropolisUpdate(LatticeType& lattice, double J, double beta, vector<double>& energy, int L, int index) {
    int pos = (int)rand() % (L * L); // pos is a random number between 0 and L^2    
    int y = pos % L;
    int x = (pos - y) / L;
    int spin_value = lattice[x][y];
    double Delta_E = evaluate_Delta_E(spin_value, lattice, J, x, y, L);

    // Metropolis step
    if (Delta_E < 0 || log((double)rand() / RAND_MAX) < -beta * Delta_E) {
        // Accept the move
        lattice[x][y] = -spin_value;
        energy[index] += Delta_E;
    }
}


int main(){
    /*
    ========================================
    = Exercise A.7: Multiple Markov Chains =
    ========================================
    */

    // Set up the rnd generator
    int iseed=3; //10112023;
    srand(static_cast<unsigned>(iseed));

    // Ising lattice
    int L = 50; // Linear dimension

    // Work with 5 different Markov Chains
    int** lattice_1 = generateRandomLattice(L);
    int** lattice_2 = generateRandomLattice(L);
    int** lattice_3 = generateRandomLattice(L);
    int** lattice_4 = generateRandomLattice(L);
    int** lattice_5 = generateRandomLattice(L);

    double J = 1.;

    int x = 0;
    int y = 0;
    int pos = 0;

    int MC_step = L*L;
    int tot_steps = pow(10, 6);
    double Delta_E = 0;
    double prob_1_2 = 0;
    double prob_2_3 = 0;
    double prob_3_4 = 0;
    double prob_4_5 = 0;
    int swap_frequency = 100/20;
    int spin_value = 0;
    double k_B = 1.;
    double T_c = 2*J/(k_B*log(1+sqrt(2)));
    double eps = 0.01; 
    vector<double> TT = {T_c-2.*eps, T_c-eps, T_c, T_c+eps, T_c+2.*eps};
    vector<double> beta = {1./(k_B*TT[0]), 1./(k_B*TT[1]), 1./(k_B*TT[2]), 
                           1./(k_B*TT[3]), 1./(k_B*TT[4])};
    
    double *Energy_observale_1 = (double*)malloc(tot_steps * sizeof(double));
    // double *Magnetization_1 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_2 = (double*)malloc(tot_steps * sizeof(double));
    // double *Magnetization_2 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_3 = (double*)malloc(tot_steps * sizeof(double));
    // double *Magnetization_3 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_4 = (double*)malloc(tot_steps * sizeof(double));
    // double *Magnetization_4 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_5 = (double*)malloc(tot_steps * sizeof(double));
    // double *Magnetization_5 = (double*)malloc(tot_steps * sizeof(double));
    
    vector<double> energy = {compute_energy(lattice_1, J, L), compute_energy(lattice_2, J, L), compute_energy(lattice_3, J, L),
                             compute_energy(lattice_4, J, L), compute_energy(lattice_5, J, L)};

    Energy_observale_1[0] = energy[0]/(L*L);
    // Magnetization_1[0] = sumMatrix(lattice_1, L)/(L*L);
    Energy_observale_2[0] = energy[1]/(L*L);
    // Magnetization_2[0] = sumMatrix(lattice_2, L)/(L*L);
    Energy_observale_3[0] = energy[2]/(L*L);
    // Magnetization_3[0] = sumMatrix(lattice_3, L)/(L*L);
    Energy_observale_4[0] = energy[3]/(L*L);
    // Magnetization_4[0] = sumMatrix(lattice_4, L)/(L*L);
    Energy_observale_5[0] = energy[4]/(L*L);
    // Magnetization_5[0] = sumMatrix(lattice_5, L)/(L*L);

    string s = "";
    // Store initial configuration
    // s = "res/MMC/T_c_2eps_lattice.txt";
    // saveMatrixToFile_int(lattice_1, L, L, s);
    // s = "res/MMC/T_c_eps_lattice.txt";
    // saveMatrixToFile_int(lattice_2, L, L, s);
    // s = "res/MMC/T_c_lattice.txt";
    // saveMatrixToFile_int(lattice_3, L, L, s);
    // s = "res/MMC/T_c_p_eps_lattice.txt";
    // saveMatrixToFile_int(lattice_4, L, L, s);
    // s = "res/MMC/T_c_p_2eps_lattice.txt";
    // saveMatrixToFile_int(lattice_5, L, L, s);
    
    
    // Create temp lattice and energy
    int** tmp_lattice = generateRandomLattice(L);

    // Keep track of swapping rates
    vector<vector<int>> swaps;


    for (int t = 0; t < tot_steps; t++){
        for (int mc = 0; mc < MC_step; mc++){
            // Select a random place on the 2D lattice
            // For the 5 chains considered:
            metropolisUpdate(lattice_1, J, beta[0], energy, L, 0);
            metropolisUpdate(lattice_2, J, beta[1], energy, L, 1);
            metropolisUpdate(lattice_3, J, beta[2], energy, L, 2);
            metropolisUpdate(lattice_4, J, beta[3], energy, L, 3);
            metropolisUpdate(lattice_5, J, beta[4], energy, L, 4);
        }

        // Store results per Metropolis sweep
        Energy_observale_1[t] = energy[0]/(L*L);
        // Magnetization_1[t] = sumMatrix(lattice_1, L)/(L*L);
        Energy_observale_2[t] = energy[1]/(L*L);
        // Magnetization_2[t] = sumMatrix(lattice_2, L)/(L*L);
        Energy_observale_3[t] = energy[2]/(L*L);
        // Magnetization_3[t] = sumMatrix(lattice_3, L)/(L*L);
        Energy_observale_4[t] = energy[3]/(L*L);
        // Magnetization_4[t] = sumMatrix(lattice_4, L)/(L*L);
        Energy_observale_5[t] = energy[4]/(L*L);
        // Magnetization_5[t] = sumMatrix(lattice_5, L)/(L*L);
    }

    
    s = "res/MMC/"+to_string(L)+"/T_c_2eps_energy.txt";
    saveResultsToFile(Energy_observale_1, tot_steps, s);
    s = "res/MMC/"+to_string(L)+"/T_c_eps_energy.txt";
    saveResultsToFile(Energy_observale_2, tot_steps, s);
    s = "res/MMC/"+to_string(L)+"/T_c_energy.txt";
    saveResultsToFile(Energy_observale_3, tot_steps, s);
    s = "res/MMC/"+to_string(L)+"/T_c_p_eps_energy.txt";
    saveResultsToFile(Energy_observale_4, tot_steps, s);
    s = "res/MMC/"+to_string(L)+"/T_c_p_2eps_energy.txt";
    saveResultsToFile(Energy_observale_5, tot_steps, s);

    fprintf(stdout, "Simulation with L = %d completed!\n", L);

    // free memory
    free(Energy_observale_1);
    // free(Magnetization_1);
    free(Energy_observale_2);
    // free(Magnetization_2);
    free(Energy_observale_3);
    // free(Magnetization_3);
    free(Energy_observale_4);
    // free(Magnetization_4);
    free(Energy_observale_5);
    // free(Magnetization_5);


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
