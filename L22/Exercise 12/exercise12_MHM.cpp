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

vector<double> implicit_func(vector<double> Z_k, vector<double> T_vals, vector<vector<double>> all_data, vector<int> corr_time_e);
double convergence(vector<double> Z_k_t, vector<double> Z_k_t_1);
vector<double> rescale(vector<double>);
vector<double> fixed_point(vector<double> T_vals, vector<vector<double>> all_data, vector<int> corr_time_e, vector<double> (*func)(vector<double>, vector<double>, vector<vector<double>>, vector<int>), vector<double> initial_guess, double tolerance, double max_iterations);

double vec_min(vector<double>);
double vec_max(vector<double>);

double compute_corr_func(vector<double> &data, int, int);

int main(){
    /*
    ==================================================
    = Exercise 1.13: Multiple Histogram Method (MHM) =
    ==================================================
    */

    int L = 50;
    double N = pow(L, 2.);

    string folder = "res/MMC/"+to_string(L)+"/";
    double T_c = 2./log(1+sqrt(2));
    double eps = 0.01;
    vector<double> T_vals = {T_c-2.*eps, T_c-eps, T_c, T_c+eps, T_c+2.*eps};
    vector<string> T_list = {"T_c_2eps_", "T_c_eps_", "T_c_", "T_c_p_eps_", "T_c_p_2eps_"};

    // string folder = "../../L10/Exercise 4/MMC_no_swap/";
    // double T_c = 2./log(1+sqrt(2));
    // vector<double> T_vals = {0.5*T_c, 0.8*T_c, T_c, 1.5*T_c, 2.*T_c};
    // vector<string> T_list = {"0.5_T_c_", "0.8_T_c_", "T_c_", "1.5_T_c_", "2_T_c_"};


    vector<string> keys = {"lattice", "magnetization", "energy", "swapping_rates"};
    string ex = ".txt";

    vector<vector<double>> all_data;
    int index = 2;
    string file_name = "";
    vector<double> average_energy;
    

    vector<int> corr_time_e;
    int t_max = pow(10, 4);

    for (int k = 0; k < T_list.size(); k++){
        string T = T_list[k];
        file_name = folder + T + keys[index] + ex;
        vector<double> data = readFromFile(file_name);
        // Passing from energy per spin to total energy
        for(int j = 0; j < data.size(); j++){
            data[j] = data[j]*N;
        }
        all_data.push_back(data);
        average_energy.push_back(sum_vector(data)/data.size());    

        // Compute correlation time
        vector<double> correlation_function_e;
        for (int t = 0; t < t_max; t++){
            correlation_function_e.push_back(compute_corr_func(data, t_max, t));
        }
        
        double normalization = correlation_function_e[0];
        for (int t = 0; t < t_max; t++){
            // rescaling by the first value
            correlation_function_e[t] /= normalization;
        }
        // saveVectorToFile(correlation_function_e, "C_corr_func_"+to_string(k)+".txt");
        corr_time_e.push_back(int(sum_vector(correlation_function_e)));
        // cout << "The correlation time for k = " << k << " is " << corr_time_e[k] << endl;      
    }
    
    // Set to zero to have convergence
    for(int l = 0; l < corr_time_e.size(); l++){
        corr_time_e[l] = 0;
    }

    saveVectorToFile(average_energy, "res/mhm/"+to_string(L)+"/direct_average.txt");

    cout << "Iteratively computing Z_k ... " << endl;
    // vector<double> Z_k(T_vals.size(), 1.0);
    vector<double> Z_k = readFromFile("res/mhm/"+to_string(L)+"/Z_k.txt");

    Z_k = fixed_point(T_vals, all_data, corr_time_e, &implicit_func, Z_k, 1e-7, 100);
    for (int k = 0; k < Z_k.size(); k++){
        cout << "k = " << k << endl;
        cout << "Z_k[k] = " << Z_k[k] << endl;
    }

    saveVectorToFile(Z_k, "res/mhm/"+to_string(L)+"/Z_k.txt");

    
    // Compute Z(β) for all βs
    vector<double> beta_j;
    vector<int> M;
    for (int j = 0; j < T_vals.size(); j++){
        beta_j.push_back(1./T_vals[j]);
        M.push_back(all_data[j].size()-corr_time_e[j]);
    }

    vector<double> Z_beta;
    vector<double> U_beta;
    vector<double> betas;

    int num = 50;
    // Select zoom if you want prediction close to the working temperatures
    bool zoom = true;
    double minimum = 0, maximum = 0;
    if(zoom){
        minimum = vec_min(beta_j);
        maximum = vec_max(beta_j);
    }
    else{        
        minimum = 0.360;
        maximum = 0.460;
    }

    for(int r = 0; r < num; r++){
        betas.push_back(minimum + r*(maximum - minimum)/(num-1));
    } 
    if(zoom){
        saveVectorToFile(betas, "res/mhm/"+to_string(L)+"/zoom/betas.txt");
    }
    else{
        saveVectorToFile(betas, "res/mhm/"+to_string(L)+"/betas.txt");
    }
    cout << "Z_beta ... " << endl;
    for(int k = 0; k < num; k++){
        double Z_b = 0;
        for(int i = 0; i < Z_k.size(); i++){
            for(int n = corr_time_e[i]; n < all_data[i].size(); n++){
                double denom = 0;
                for(int j = 0; j < Z_k.size(); j++){
                    denom += M[j]*exp((betas[k]-beta_j[j])*all_data[i][n])/Z_k[j];
                }
                Z_b += 1./denom;
            }
        }
        Z_beta.push_back(Z_b);
    }
    if(zoom){
        saveVectorToFile(Z_beta, "res/mhm/"+to_string(L)+"/zoom/Z_beta.txt");
    }
    else{
        saveVectorToFile(Z_beta, "res/mhm/"+to_string(L)+"/Z_beta.txt");
    }
    // Compute U(β) for all βs
    cout << "U_beta ... " << endl;
    for(int k = 0; k < num; k++){
        double U_b = 0;
        for(int i = 0; i < Z_k.size(); i++){
            for(int n = corr_time_e[i]; n < all_data[i].size(); n++){
                double denom = 0;
                for(int j = 0; j < Z_k.size(); j++){
                    denom += M[j]*exp((betas[k]-beta_j[j])*all_data[i][n])/Z_k[j];
                }
                U_b += all_data[i][n]/denom;
            }
        }
        U_beta.push_back(U_b/Z_beta[k]);
    }
    if(zoom){    
        saveVectorToFile(U_beta, "res/mhm/"+to_string(L)+"/zoom/U_beta.txt");
    }
    else{
        saveVectorToFile(U_beta, "res/mhm/"+to_string(L)+"/U_beta.txt");
    }
    // Compute C(β) for all βs
    // Fist compute <E^2>_β
    vector<double> spec_heat;
    cout << "C_beta ... " << endl;
    for(int k = 0; k < num; k++){
        double E_2_b = 0;
        for(int i = 0; i < Z_k.size(); i++){
            for(int n = corr_time_e[i]; n < all_data[i].size(); n++){
                double denom = 0;
                for(int j = 0; j < Z_k.size(); j++){
                    denom += M[j]*exp((betas[k]-beta_j[j])*all_data[i][n])/Z_k[j];
                }
                E_2_b += pow(all_data[i][n], 2.)/denom;
            }
        }
        E_2_b /= Z_beta[k];
        // (<E>_β)^2 is given from the previous computation
        spec_heat.push_back((E_2_b-pow(U_beta[k], 2.))/(N*pow(betas[k], 2.)));
    }
    if (zoom){ 
        saveVectorToFile(spec_heat, "res/mhm/"+to_string(L)+"/zoom/C_beta.txt");
    }
    else{
        saveVectorToFile(spec_heat, "res/mhm/"+to_string(L)+"/C_beta.txt");
    }

    return 0;
}


double sum_vector(vector<double> &x){
    double res = 0.;
    for(auto val : x){
        res += val;
    }
    return res;
}

vector<double> implicit_func(vector<double> Z_k, vector<double> T_vals, vector<vector<double>> all_data, vector<int> corr_time_e){
    vector<double> out_Z_k(T_vals.size(), 0.); 
    vector<double> beta_j;
    vector<int> M;
    for (int j = 0; j < T_vals.size(); j++){
        // if (T_vals[j] != T_vals[k]){
        beta_j.push_back(1./T_vals[j]);
        M.push_back(all_data[j].size()-corr_time_e[j]);
        // }
    }
    for (int k = 0; k < T_vals.size(); k++){
        double beta_k = 1./T_vals[k];
        
        // Start computation
        double new_Z_k = 0;
        for (int i = 0; i < beta_j.size(); i++){
            for(int n = corr_time_e[i]; n < all_data[i].size(); n++){
                double denom = 0;
                for(int j = 0; j < M.size(); j++){
                    denom += M[j]*exp((beta_k-beta_j[j])*all_data[i][n])/Z_k[j];
                }
                new_Z_k += 1./denom;
            }
        }

        out_Z_k[k] = new_Z_k;

    }
    
    return out_Z_k;
}


double convergence(vector<double> Z_k_t, vector<double> Z_k_t_1){
    vector<double> delta;
    for(int k = 0; k < Z_k_t.size(); k++){
        delta.push_back(pow((Z_k_t[k]-Z_k_t_1[k])/Z_k_t[k], 2.0));
    }
    return sum_vector(delta);
}


vector<double> rescale(vector<double> x){
    double minimum = vec_min(x);
    double maximum = vec_max(x);
    double A = 1./sqrt(minimum*maximum);
    vector<double> x_new;
    for(auto val : x){
        x_new.push_back(val*A);
    }
    return x_new; 
}


vector<double> fixed_point(vector<double> T_vals, vector<vector<double>> all_data, vector<int> corr_time_e, vector<double> (*func)(vector<double>, vector<double>, vector<vector<double>>, vector<int>), vector<double> initial_guess, double tolerance, double max_iterations){
    vector<double> x = initial_guess;
    vector<double> x_new;
    for(int i = 0; i < max_iterations; i++){
        x_new = func(x, T_vals, all_data, corr_time_e);
        // Rescale Z_k to avoid underflows
        if (convergence(x_new, x) < tolerance*tolerance){
            cout << "Convergence at iteration " << i << "\n" << endl;
            return x_new;
        }
        // x_new = rescale(x_new);
        x = x_new; 
    }
    cout << "Fixed-point iteration did not converge within the specified max number of iterations.\n" << endl;
    return x;
}

double vec_min(vector<double> v){
    double res = v[0];
    for(auto a : v){
        if (a < res){
            res = a;
        }
    }
    return res;
}

double vec_max(vector<double> v){
    double res = v[0];
    for(auto a : v){
        if (a > res){
            res = a;
        }
    }
    return res;
}

double compute_corr_func(vector<double> &data, int t_max, int t){
    double res = 0;
    for(int o = 0; o < t_max-t; o++){
        res += data[o]*data[t+o];
    }
    res /= (t_max-t);
    
    vector<double> data_1(data.begin(), data.begin()+(t_max-t));
    vector<double> data_2(data.begin()+t, data.begin()+t_max);
    res += -(sum_vector(data_1) / (t_max - t)) * (sum_vector(data_2) / (t_max - t));

    return res;
}