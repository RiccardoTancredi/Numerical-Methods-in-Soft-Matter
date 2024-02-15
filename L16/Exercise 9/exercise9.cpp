#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

/* our custom libraries */
#include "save_data.h" // Include the custom header file

double Euler_integrator(double dt, double x_t, double p_t, string type);
double Symplectic_integrator(double dt, double x_t, double p_t, string type);
double Exact_solution(double dt, double x_t, double p_t, string type);
double sumMatrix(int **matrix, int L);


int main(){
    /*
    ======================================
    = Exercise 1.10: Integration schemes =
    ======================================
    */

    // Set up the rnd generator
    int iseed=3; //12012024;
    srand(static_cast<unsigned>(iseed));

    // integration time
    vector<double> delta_t = {pow(10, -3.), pow(10, -2.)}; 
    double x_0 = 1., p_0 = 0.;

    vector<double> x_Euler;
    vector<double> p_Euler;  
    vector<double> x_Symp;
    vector<double> p_Symp;  
    vector<double> x_Exact;
    vector<double> p_Exact;  

    int T = 10;
    string file_path = "res/";

    for(int i = 0; i < delta_t.size(); i++){
        cout << "Buffering dt = " << delta_t[i] << endl;
        string result = to_string(delta_t[i]);
        size_t found = result.find_last_not_of('0');
        if (found != std::string::npos) result.erase(found + 1, std::string::npos);
        
        // Set initial conditions
        x_Euler.push_back(x_0); p_Euler.push_back(p_0);
        x_Symp.push_back(x_0); p_Symp.push_back(p_0);
        x_Exact.push_back(x_0); p_Exact.push_back(p_0);
        for(int iter = 0; iter < T / delta_t[i]; iter++){

            // First integrator
            x_Euler.push_back(Euler_integrator(delta_t[i], x_Euler[iter], p_Euler[iter], "x"));
            p_Euler.push_back(Euler_integrator(delta_t[i], x_Euler[iter], p_Euler[iter], "p"));

            // Second integrator
            x_Symp.push_back(Symplectic_integrator(delta_t[i], x_Symp[iter], p_Symp[iter], "x"));
            p_Symp.push_back(Symplectic_integrator(delta_t[i], x_Symp[iter], p_Symp[iter], "p"));
            
            // Exact solution
            x_Exact.push_back(Exact_solution((iter+1)*delta_t[i], x_0, p_0, "x"));
            p_Exact.push_back(Exact_solution((iter+1)*delta_t[i], x_0, p_0, "p"));  
        }

        cout << "Saving results..." << endl;

        // Saving results...
        saveVectorToFile(x_Euler, file_path + result + "/Euler_x.dat");
        saveVectorToFile(p_Euler, file_path + result + "/Euler_p.dat");

        saveVectorToFile(x_Symp, file_path + result + "/Symp_x.dat");
        saveVectorToFile(p_Symp, file_path + result + "/Symp_p.dat");
        
        saveVectorToFile(x_Exact, file_path + result + "/Exact_x.dat");
        saveVectorToFile(p_Exact, file_path + result + "/Exact_p.dat");

        // Clear lists
        x_Euler.clear(); p_Euler.clear(); 
        x_Symp.clear(); p_Symp.clear();
        x_Exact.clear(); p_Exact.clear();
    }


    return 0;
}

double Euler_integrator(double dt, double x_t, double p_t, string type){
    if (type == "x"){
        return x_t + p_t * dt;
    }
    else if(type == "p"){
        return p_t - x_t * dt;
    }
    else{
        cout << "Select a variable to integrate: x or p!" << endl;
        exit(1);
    }
}
double Symplectic_integrator(double dt, double x_t, double p_t, string type){
    if (type == "x"){
        return x_t * (1.-dt*dt) + p_t * dt;
    }
    else if(type == "p"){
        return p_t - x_t * dt;
    }
    else{
        cout << "Select a variable to integrate: x or p!" << endl;
        exit(1);
    }
}
double Exact_solution(double dt, double x_0, double p_0, string type){
    if (type == "x"){
        return cos(dt)* x_0 + sin(dt) * p_0;
    }
    else if(type == "p"){
        return -sin(dt) * x_0 + cos(dt) * p_0;
    }
    else{
        cout << "Select a variable to integrate: x or p!" << endl;
        exit(1);
    }
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
