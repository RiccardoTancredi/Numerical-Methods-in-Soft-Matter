#define N_MOVES 1
#define PI 3.141592653589793

typedef struct
{
    double x, y, z;     // position

    double vx, vy, vz;  // velocity

    double fx, fy, fz;  // force

    double ofx, ofy, ofz;  // force

    double ox, oy, oz;	//old colloid coordinates
    
    double ovx, ovy, ovz;	//old colloid velocities
    
    double tx, ty, tz;  //temp colloid coordinates

    double Cx, Cy, Cz;  // C variables of the Underdamped Langevin eq.

} Particle;

Particle *parts = NULL; //It is always a good practice to assign a NULL value to a pointer variable in case you do not have an exact address to be assigned.

//defining the System of patchy IPCs properties

typedef struct  
{
    int NPart;
    int step;
    int NSteps, NPrint, Nrestart;      //number of IPCs, number of MC sweeps, frequency of configurations
    int restart;
    long seed;

    int overlap;
    long int tries[N_MOVES];		//array with trial number for different moves
    long int accepted[N_MOVES];		//acceptance for different moves

    double disp_max;		       //maximum displacement for colloids
	
    double box_x, box_y, box_z;
    double T;	    // temperature
    double gamma;   // friction coefficient	
    double K;       // spring coefficient	
    double energy;
    
    double sigma;
    double m;   // mass
    double sigma_cut;
    double eps;
    int accept;

    int model;			            // switch for models
    int integrator;			        // switch for integrator
    int exercise;			        // switch exercise
    double spring_2_x;			    // position of the second spring: x
    double spring_2_y;			    // position of the second spring: y
    char start_file[10]; 
    char velocity_file[20]; 
    char kind[10];
    
    int iter;                       // iteration: select the folder where to store the simulation data


} System;

System mySys;

