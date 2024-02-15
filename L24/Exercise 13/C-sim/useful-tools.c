#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "array_alloc.h"
#include "useful-tools.h"


// WARNING - inline function was removed from .h because it did not compile well

/* https://gcc.gnu.org/onlinedocs/gcc/Inline.html :
   By declaring a function inline, you can direct GCC to make calls to that function faster. One way GCC can achieve this is to integrate that function's code into the code for its callers. This makes execution faster by eliminating the function-call overhead */
inline int Periodic_int(int a, int A){
  return (a % A + A) % A;
}

inline double Periodic(double x, double X){
  while(x<0)  x+=X;
  while(x>=X) x-=X;
  return x;
}

inline double Shortest(double v, double X_2, double X){
  if(v > X_2)  return v-X;
  if(v < -X_2) return v+X;
  return v;
}

// Read parameters from a file -- accelarate compilation

void readParameters(const char *filename, int *nstep_save, int *N, int *Lp, double *box, double *T,
                    double *v_trap_ini, int *Nv, double *v_per_decade, double *k_trap, double *k_pol,
                    double *R, double *eps, double *R0, double *eps0, double *dt, double *tt, 
                    double *Lambda, double *f_active, int *model)
{
    FILE *_myfile = fopen(filename, "r");
    
    char line[100];
    int i=0;
    char *token;
    char *_key=NULL;
    char *_value = NULL;

    while(fgets(line, sizeof(line), _myfile) != NULL){

      token = strtok(line, " ");
      i = 0;
      while (token != NULL) {
          if (i == 0) _key = token;
	  if (i == 2) _value = token;
          token = strtok(NULL, " ");
          i++;
      }     

      /*
      saving configs if > 0
      N = # of small particles
      L = length of each polymer
      box[0] = box size, x
      box[1] = box size, y
      T = temperature (kB=1)
      v_trap_ini = smallest velocity of the trap: NOT USED
      Nv = number of velocities: NOT USED
      v_per_decade = velocities per decade: NOT USED
      k_trap = stiffness of the trap 
      k_pol = stiffness of the polymer bonds
      R = "radius" of each particle, but if(R>0.25) -> STOP PROGRAMM
      eps = repulsive energy of particles
      R0 = "radius" of the probe
      eps0 = repulsive energy of the probe
      dt = integration time step
      tt = total time of the simulation
    */

      if(strcmp("nstep_save",_key) == 0){ *nstep_save = atoi(_value); continue; }
      if(strcmp("N",_key) == 0){ *N = atoi(_value); continue; }
      if(strcmp("Lp",_key) == 0){ *Lp = atoi(_value); continue; }
      if(strcmp("box_x",_key) == 0){ box[0] = atof(_value); continue; } 
      if(strcmp("box_y",_key) == 0){ box[1] = atof(_value); continue; }
      // if(strcmp("box_z",_key) == 0){ box[2] = atof(_value); continue; }
      if(strcmp("T",_key) == 0){ *T = atof(_value); continue; }
      if(strcmp("v_trap_ini",_key) == 0){ *v_trap_ini = atof(_value); continue; }
      if(strcmp("Nv",_key) == 0){ *Nv = atoi(_value); continue; } 
      if(strcmp("v_per_decade",_key) == 0){ *v_per_decade = atof(_value); continue; }
      if(strcmp("k_trap",_key) == 0){ *k_trap = atof(_value); continue; }
      if(strcmp("k_pol",_key) == 0){ *k_pol = atof(_value); continue; }
      if(strcmp("R",_key) == 0){ *R = atof(_value); continue; }
      if(strcmp("eps",_key) == 0){ *eps = atof(_value); continue; }    
      if(strcmp("R0",_key) == 0){ *R0 = atof(_value); continue; }
      if(strcmp("eps0",_key) == 0){ *eps0 = atof(_value); continue; }
      if(strcmp("dt",_key) == 0){ *dt = atof(_value); continue; }
      if(strcmp("tt",_key) == 0){ *tt = atoi(_value); continue; }         
      if(strcmp("Lambda",_key) == 0){ *Lambda = atof(_value); continue; }         
      if(strcmp("f_active",_key) == 0){ *f_active = atof(_value); continue; }     
      if(strcmp("model",_key) == 0){ *model = atoi(_value); continue; }    
    }
    
    fclose(_myfile);

}

// void readParameters(const char *filename, int *nstep_save, int *N, int *Lp, double *box, double *T,
//                     double *v_trap_ini, int *Nv, double *v_per_decade, double *k_trap, double *k_pol,
//                     double *R, double *eps, double *R0, double *eps0, double *dt, double *tt) {
    
//     /*
//       saving configs if > 0
//       N = nr of small particles
//       L = length of each polymer
//       box[0] = box size, x
//       box[1] = box size, y
//       T = temperature (kB=1)
//       v_trap_ini = smallest velocity of the trap: NOT USED
//       Nv = number of velocities: NOT USED
//       v_per_decade = velocities per decade: NOT USED
//       k_trap = stiffness of the trap 
//       k_pol = stiffness of the polymer bonds
//       R = "radius" of each particle, but if(R>0.25) -> STOP PROGRAMM
//       eps = repulsive energy of particles
//       R0 = "radius" of the probe
//       eps0 = repulsive energy of the probe
//       dt = integration time step
//       tt = total time of the simulation
//     */

//     FILE *file = fopen(filename, "r");
//     if (file == NULL) {
//         fprintf(stderr, "Error opening file: %s\n", filename);
//         exit(EXIT_FAILURE);
//     }

//     // Use fscanf to read parameters from the file
//     fscanf(file, "%d %d %d %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
//            nstep_save, N, Lp, &box[0], &box[1], T, v_trap_ini, Nv, v_per_decade,
//            k_trap, k_pol, R, eps, R0, eps0, dt, tt);

//     fclose(file);
// }