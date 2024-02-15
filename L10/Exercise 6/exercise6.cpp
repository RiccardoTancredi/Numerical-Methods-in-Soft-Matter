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
// double evaluate_Delta_E(int spin_value, int **lattice, double J, int x, int y, int L);
double sumMatrix(int **matrix, int L);

int main(){
    /*
    =================================
    = Exercise A.7: Wolff algorithm =
    =================================
    */

    // Set up the rnd generator
    int iseed=3; //10112023;
    srand(static_cast<unsigned>(iseed));

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
            lattice[i][j] = (double)rand() / (RAND_MAX) < 0.5 ? -1:1;
        }
    }

    double J = 1.;

    int x = 0;
    int y = 0;
    int pos = 0;

    // int MC_step = L*L;
    int tot_steps = pow(10, 3);
    // double Delta_E = 0;
    int spin_value = 0;
    double k_B = 1.;
    double T_c = 2*J/(k_B*log(1+sqrt(2))); 
    double T = 2*T_c;
    double beta = 1./(k_B*T);
    double P_add = 1 - exp(-2*beta*J);

    double *Energy_observale = (double*)malloc(tot_steps * sizeof(double));
    double *Magnetization = (double*)malloc(tot_steps * sizeof(double));
    double *Cluster_size = (double*)malloc(tot_steps * sizeof(double));
    
    double energy = compute_energy(lattice, J, L);

    Energy_observale[0] = energy/(L*L);
    Magnetization[0] = sumMatrix(lattice, L)/(L*L);

    cout << "The P_add probability is = " << P_add << endl;
    vector<int> neigh_pos(4, 0);

    set<int> neighbours;
    set<int> new_neighbours;
    set<int> cluster;

    string s = "res/0_2_T_c_lattice.txt";
    saveMatrixToFile_int(lattice, L, L, s);
    int spin_pos = 0;
    int prev_length = 0;

    for (int t = 1; t < tot_steps; t++){
        // Select a random place on the 2D lattice
        pos = (int)rand() % (L*L);   // pos is a random number between 0 and L^2    
        if (t % 100 == 0){
            cout << "The iteration is " << t << endl;
        }
        int y = pos % L;
        int x = (pos-y)/L;
        spin_value = lattice[x][y];
        neighbours.insert(pos);
        bool it = true;
        while(it){
            prev_length = cluster.size();
            for(auto spin_pos : neighbours){
                // Collect the 4 neighbours of spin_pos:
                neigh_pos[0] = (spin_pos-L < 0) ? spin_pos + L*(L-1): spin_pos-L;
                neigh_pos[1] = (spin_pos % L == 0) ? spin_pos + (L-1): spin_pos-1;
                neigh_pos[2] = (spin_pos % L == L-1) ? spin_pos - (L-1): spin_pos+1;
                neigh_pos[3] = (spin_pos % (L*(L-1)) < L && spin_pos > L) ? spin_pos - L*(L-1): spin_pos+L;
                
                // Start to create the cluster: look at neighbours                
                for (auto neigh : neigh_pos){
                    int y = neigh % L;
                    int x = (neigh-y)/L;
                    auto f = cluster.find(neigh);
                    if((lattice[x][y] == spin_value) && ((double)rand() / (RAND_MAX) < P_add && f == cluster.end())){
                        // Same sign of each element in the cluster: 
                        // add to cluster with probability P_add
                        new_neighbours.insert(neigh);
                    }
                }
            }
            for (int n : neighbours){
                cluster.insert(n);
            }
            neighbours.clear();
            for(int nn : new_neighbours){
                neighbours.insert(nn);
                cluster.insert(nn);
            }
            new_neighbours.clear();
            // Add new_neighbours in neighbours
            // Check if the size of the vector has increased, otherwise stop
            it = (prev_length == cluster.size()) ? false : true;
        }

        // Flip all the elements in the cluster:
        for(int el : cluster){
            int y = el % L;
            int x = (el-y)/L;
            lattice[x][y] *= -1;
        }
        
        Cluster_size[t] = cluster.size();
        energy = compute_energy(lattice, J, L);
        Energy_observale[t] = energy/(L*L);
        Magnetization[t] = sumMatrix(lattice, L)/(L*L);

        if(t % 100 == 0){
            s = "res/" + to_string(t) + "_2_T_c_lattice.txt";
            saveMatrixToFile_int(lattice, L, L, s);
        }

        // Erase all neighbours sets
        cluster.clear();
        neighbours.clear();
        new_neighbours.clear();
    }

    s = "res/" + to_string(tot_steps) + "_2_T_c_lattice.txt";
    saveMatrixToFile_int(lattice, L, L, s);
    
    s = "res/2_T_c_magnetization.txt";
    saveResultsToFile(Magnetization, tot_steps, s);

    s = "res/2_T_c_energy.txt";
    saveResultsToFile(Energy_observale, tot_steps, s);

    s = "res/2_T_c_cluster_size.txt";
    saveResultsToFile(Cluster_size, tot_steps, s);

    fprintf(stdout, "Simulation with L = %d completed!\n", L);

    // free memory
    free(lattice);
    free(Energy_observale);
    free(Magnetization);
    free(Cluster_size);

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


// double evaluate_Delta_E(int spin_value, int **lattice, double J, int x, int y, int L){
//     int neighbours_x[4] = {(x-1+L) % L, x, x, (x+1) % L};
//     int neighbours_y[4] = {y, (y-1+L) % L, (y+1) % L, y};
//     int sum_spins = 0;
//     for (int index = 0; index < 4; index++){
//         sum_spins += lattice[neighbours_x[index]][neighbours_y[index]];
//     }
//     return(2*spin_value*J*sum_spins);
// }

double sumMatrix(int **matrix, int L) {
    double sum_elements = 0.;
    for(int i = 0; i < L; i++) {
        for(int j = 0; j < L; j++) {
            sum_elements += matrix[i][j];
        }
    }
    return sum_elements;
}
