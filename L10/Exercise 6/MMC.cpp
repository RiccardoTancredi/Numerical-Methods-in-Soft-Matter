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
    vector<double> TT = {T_c/2., 0.8*T_c, T_c, 1.5*T_c, 2.*T_c};
    vector<double> beta = {1./(k_B*TT[0]), 1./(k_B*TT[1]), 1./(k_B*TT[2]), 
                           1./(k_B*TT[3]), 1./(k_B*TT[4])};
    
    double *Energy_observale_1 = (double*)malloc(tot_steps * sizeof(double));
    double *Magnetization_1 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_2 = (double*)malloc(tot_steps * sizeof(double));
    double *Magnetization_2 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_3 = (double*)malloc(tot_steps * sizeof(double));
    double *Magnetization_3 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_4 = (double*)malloc(tot_steps * sizeof(double));
    double *Magnetization_4 = (double*)malloc(tot_steps * sizeof(double));
    double *Energy_observale_5 = (double*)malloc(tot_steps * sizeof(double));
    double *Magnetization_5 = (double*)malloc(tot_steps * sizeof(double));
    
    vector<double> energy = {compute_energy(lattice_1, J, L), compute_energy(lattice_2, J, L), compute_energy(lattice_3, J, L),
                             compute_energy(lattice_4, J, L), compute_energy(lattice_5, J, L)};

    Energy_observale_1[0] = energy[0]/(L*L);
    Magnetization_1[0] = sumMatrix(lattice_1, L)/(L*L);
    Energy_observale_2[0] = energy[1]/(L*L);
    Magnetization_2[0] = sumMatrix(lattice_2, L)/(L*L);
    Energy_observale_3[0] = energy[2]/(L*L);
    Magnetization_3[0] = sumMatrix(lattice_3, L)/(L*L);
    Energy_observale_4[0] = energy[3]/(L*L);
    Magnetization_4[0] = sumMatrix(lattice_4, L)/(L*L);
    Energy_observale_5[0] = energy[4]/(L*L);
    Magnetization_5[0] = sumMatrix(lattice_5, L)/(L*L);


    // Store initial configuration
    string s = "MMC/0.5_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_1, L, L, s);
    s = "MMC/0.8_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_2, L, L, s);
    s = "MMC/T_c_lattice.txt";
    saveMatrixToFile_int(lattice_3, L, L, s);
    s = "MMC/1.5_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_4, L, L, s);
    s = "MMC/2_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_5, L, L, s);
    
    
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

        // Add swapping rate possibility between chains
        if ((t+1) % swap_frequency == 0){
            int firstNumber = (int)rand() % 5 + 1;
            int secondNumber;
            do {
                secondNumber = (int)rand() % 5 + 1;
            } while (secondNumber == firstNumber || abs(firstNumber-secondNumber) > 1);  

            prob_1_2 = (beta[1]-beta[0])*(energy[1]-energy[0]);
            prob_2_3 = (beta[2]-beta[1])*(energy[2]-energy[1]);
            prob_3_4 = (beta[3]-beta[2])*(energy[3]-energy[2]);
            prob_4_5 = (beta[4]-beta[3])*(energy[4]-energy[3]);

            if ((firstNumber == 1 && secondNumber == 2) || (firstNumber == 2 && secondNumber == 1)){
                if (prob_1_2 < 0 || log((double)rand() / RAND_MAX) < prob_1_2) {
                    // Accept the move
                    tmp_lattice = lattice_1;
                    lattice_1 = lattice_2;
                    lattice_2 = tmp_lattice;
                    swap(energy[0], energy[1]);
                    swaps.push_back({firstNumber, secondNumber});
                }
            }
            else if ((firstNumber == 2 && secondNumber == 3) || (firstNumber == 3 && secondNumber == 2)){
                if (prob_2_3 < 0 || log((double)rand() / RAND_MAX) < prob_2_3) {
                    // Accept the move
                    tmp_lattice = lattice_2;
                    lattice_2 = lattice_3;
                    lattice_3 = tmp_lattice;
                    swap(energy[1], energy[2]);
                    swaps.push_back({firstNumber, secondNumber});
                }
            }
            else if ((firstNumber == 3 && secondNumber == 4) || (firstNumber == 4 && secondNumber == 3)){
                if (prob_3_4 < 0 || log((double)rand() / RAND_MAX) < prob_3_4) {
                    // Accept the move
                    tmp_lattice = lattice_3;
                    lattice_3 = lattice_4;
                    lattice_4 = tmp_lattice;
                    swap(energy[2], energy[3]);
                    swaps.push_back({firstNumber, secondNumber});
                }
            }
            else if ((firstNumber == 4 && secondNumber == 5) || (firstNumber == 5 && secondNumber == 4)){
                if (prob_4_5 < 0 || log((double)rand() / RAND_MAX) < prob_4_5) {
                    // Accept the move
                    tmp_lattice = lattice_4;
                    lattice_4 = lattice_5;
                    lattice_5 = tmp_lattice;
                    swap(energy[3], energy[4]);
                    swaps.push_back({firstNumber, secondNumber});
                }
            }
        }

        // Store results per Metropolis sweep
        Energy_observale_1[t] = energy[0]/(L*L);
        Magnetization_1[t] = sumMatrix(lattice_1, L)/(L*L);
        Energy_observale_2[t] = energy[1]/(L*L);
        Magnetization_2[t] = sumMatrix(lattice_2, L)/(L*L);
        Energy_observale_3[t] = energy[2]/(L*L);
        Magnetization_3[t] = sumMatrix(lattice_3, L)/(L*L);
        Energy_observale_4[t] = energy[3]/(L*L);
        Magnetization_4[t] = sumMatrix(lattice_4, L)/(L*L);
        Energy_observale_5[t] = energy[4]/(L*L);
        Magnetization_5[t] = sumMatrix(lattice_5, L)/(L*L);
    }

    // Store final results
    s = "MMC/" + to_string(tot_steps) + "_0.5_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_1, L, L, s);
    s = "MMC/" + to_string(tot_steps) + "_0.8_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_2, L, L, s);
    s = "MMC/" + to_string(tot_steps) + "_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_3, L, L, s);
    s = "MMC/" + to_string(tot_steps) + "_1.5_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_4, L, L, s);
    s = "MMC/" + to_string(tot_steps) + "_2_T_c_lattice.txt";
    saveMatrixToFile_int(lattice_5, L, L, s);
    

    s = "MMC/0.5_T_c_magnetization.txt";
    saveResultsToFile(Magnetization_1, tot_steps, s);
    s = "MMC/0.8_T_c_magnetization.txt";
    saveResultsToFile(Magnetization_2, tot_steps, s);
    s = "MMC/T_c_magnetization.txt";
    saveResultsToFile(Magnetization_3, tot_steps, s);
    s = "MMC/1.5_T_c_magnetization.txt";
    saveResultsToFile(Magnetization_4, tot_steps, s);
    s = "MMC/2_T_c_magnetization.txt";
    saveResultsToFile(Magnetization_5, tot_steps, s);

    s = "MMC/0.5_T_c_energy.txt";
    saveResultsToFile(Energy_observale_1, tot_steps, s);
    s = "MMC/0.8_T_c_energy.txt";
    saveResultsToFile(Energy_observale_2, tot_steps, s);
    s = "MMC/T_c_energy.txt";
    saveResultsToFile(Energy_observale_3, tot_steps, s);
    s = "MMC/1.5_T_c_energy.txt";
    saveResultsToFile(Energy_observale_4, tot_steps, s);
    s = "MMC/2_T_c_energy.txt";
    saveResultsToFile(Energy_observale_5, tot_steps, s);

    s = "MMC/swapping_rates.txt";
    ofstream outputFile(s);
    if (outputFile.is_open()) {
        for (const auto& row : swaps) {
            for (const auto& value : row) {
                outputFile << value << ' ';
            }
            outputFile << '\n';
        }
    outputFile.close();
    }

    fprintf(stdout, "Simulation with L = %d completed!\n", L);

    // free memory
    free(Energy_observale_1);
    free(Magnetization_1);
    free(Energy_observale_2);
    free(Magnetization_2);
    free(Energy_observale_3);
    free(Magnetization_3);
    free(Energy_observale_4);
    free(Magnetization_4);
    free(Energy_observale_5);
    free(Magnetization_5);


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
