#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

/* our custom libraries */
#include "save_data.h" // Include the custom header file

void Velocity_Verlet(double &x_t, double &p_t, double f_t, double omega, double dt);
double Beeman_integrator(double dt, double x_t, double p_t, double f_t, double f_t_m_dt, double f_t_p_dt, string type);
double Exact_solution(double dt, double x_t, double p_t, double omega, string type);
double sumMatrix(int **matrix, int L);


int main(){
    /*
    =================================================================
    = Exercise 1.10: Integration schemes - Velocity Verlet & others =
    =================================================================
    */

    // Set up the rnd generator
    int iseed=3; //12012024;
    srand(static_cast<unsigned>(iseed));

    double x_0 = 1., p_0 = 0.;
    double omega = 1.;
    // integration time
    double dt = pow(10, -3.);

    vector<double> x_V_V;
    vector<double> p_V_V;  
    vector<double> x_Beeman;
    vector<double> p_Beeman;  
    vector<double> x_Exact;
    vector<double> p_Exact;  

    int T = 10;
    string file_path = "res/";

    string result = to_string(dt);
    size_t found = result.find_last_not_of('0');
    if (found != std::string::npos) result.erase(found + 1, std::string::npos);
    
    // Set initial conditions
    x_V_V.push_back(x_0); p_V_V.push_back(p_0);
    x_Beeman.push_back(x_0); p_Beeman.push_back(p_0);
    x_Exact.push_back(x_0); p_Exact.push_back(p_0);
        
    // Need another point for Beeman integrator
    double x_t_m_dt = Exact_solution(-dt, x_0, p_0, omega, "x");
    vector<double> f;
    f.push_back(- x_t_m_dt * omega * omega); f.push_back(- x_0 * omega * omega); 
    
    for(int iter = 0; iter < T / dt; iter++){

        double x_t = x_V_V[iter]; double p_t = p_V_V[iter];
        double f_t = - x_t * omega * omega;
        // Velocity Verlet integrator
        Velocity_Verlet(x_t, p_t, f_t, omega, dt);
        x_V_V.push_back(x_t); p_V_V.push_back(p_t);
        
        // Beeman integrator
        x_Beeman.push_back(Beeman_integrator(dt, x_Beeman[iter], p_Beeman[iter], f[iter+1], f[iter], 0, "x"));
        // update forces
        f.push_back(- x_Beeman[iter+1] * omega * omega);
        p_Beeman.push_back(Beeman_integrator(dt, x_Beeman[iter], p_Beeman[iter], f[iter+1], f[iter], f[iter+2], "p"));
        
        // Exact solution
        x_Exact.push_back(Exact_solution((iter+1)*dt, x_0, p_0, omega, "x"));
        p_Exact.push_back(Exact_solution((iter+1)*dt, x_0, p_0, omega, "p"));  
    }

    cout << "Saving results..." << endl;

    // Saving results...
    saveVectorToFile(x_V_V, file_path + result + "/Velocity_Verlet_x.dat");
    saveVectorToFile(p_V_V, file_path + result + "/Velocity_Verlet_p.dat");

    saveVectorToFile(x_Beeman, file_path + result + "/Beeman_x.dat");
    saveVectorToFile(p_Beeman, file_path + result + "/Beeman_p.dat");
    
    saveVectorToFile(x_Exact, file_path + result + "/Exact_omega_x.dat");
    saveVectorToFile(p_Exact, file_path + result + "/Exact_omega_p.dat");

    // Clear lists
    x_V_V.clear(); p_V_V.clear(); 
    x_Beeman.clear(); p_Beeman.clear();
    x_Exact.clear(); p_Exact.clear();


    return 0;
}

double Exact_solution(double dt, double x_0, double p_0, double omega, string type){
    if (type == "x"){
        return cos(dt*omega)* x_0 + sin(dt*omega) * p_0;
    }
    else if(type == "p"){
        return -sin(dt*omega) * omega * x_0 + cos(dt*omega) * p_0 * omega;
    }
    else{
        cout << "Select a variable to integrate: x or p!" << endl;
        exit(1);
    }
}


void Velocity_Verlet(double &x_t, double &p_t, double f_t, double omega, double dt){
    
    // momentum translation by an amount (half a kick)
    p_t += dt*f_t/2.;

    // position translation by an amount (drift)
    x_t += p_t * dt;

    // forces update: for this specific model, F = -rω^2 (harmonic oscillator)
    f_t = - x_t * omega * omega;

    // momentum translation by an amount (half a kick)
    p_t += dt*f_t/2.;

}

double Beeman_integrator(double dt, double x_t, double p_t, double f_t, double f_t_m_dt, double f_t_p_dt, string type){

    if (type == "x"){
        return x_t + p_t * dt + (4.*f_t - f_t_m_dt)*dt*dt/6.;
    }
    else if(type == "p"){
        return p_t + (2.*f_t_p_dt + 5.*f_t - f_t_m_dt)*dt/6.;
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
