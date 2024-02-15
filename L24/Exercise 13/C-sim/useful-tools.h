int Periodic_int(int a, int A);

double Periodic(double x, double X);

double Shortest(double v, double X_2, double X);

// Function to read parameters from a file
void readParameters(const char *filename, int *nstep_save, int *N, int *Lp, double *box, double *T,
                    double *v_trap_ini, int *Nv, double *v_per_decade, double *k_trap, double *k_pol,
                    double *R, double *eps, double *R0, double *eps0, double *dt, double *tt, 
                    double *Lambda, double *f_active, int *model);

