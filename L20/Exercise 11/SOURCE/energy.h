double distance_updated(Particle p, Particle q){
    // compute new updated distance
    double dxu = p.x - q.x;
    double dyu = p.y - q.y;
    double dzu = p.z - q.z;
    double dx = dxu - round(dxu/mySys.box_x)*mySys.box_x;
    double dy = dyu - round(dyu/mySys.box_y)*mySys.box_y;
    double dz = dzu - round(dzu/mySys.box_z)*mySys.box_z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double distance_rw(Particle p, Particle q){
    // compute distance
    double dxu = p.ox - q.ox;
    double dyu = p.oy - q.oy;
    double dzu = p.oz - q.oz;
    double dx = dxu - round(dxu/mySys.box_x)*mySys.box_x;
    double dy = dyu - round(dyu/mySys.box_y)*mySys.box_y;
    double dz = dzu - round(dzu/mySys.box_z)*mySys.box_z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}


void compute_diff_distance_vector(double *dxu, double *dyu, double *dzu){
    // compute difference
    *dxu = *dxu - round(*dxu/mySys.box_x)*mySys.box_x;
    *dyu = *dyu - round(*dyu/mySys.box_y)*mySys.box_y;
    *dzu = *dzu - round(*dzu/mySys.box_z)*mySys.box_z;
}



double LJ_potential_per_particle_correction(){
    return -4.*mySys.eps*(pow(mySys.sigma/mySys.sigma_cut, 12.)-pow(mySys.sigma/mySys.sigma_cut, 6.));
}


double LJ_potential_per_particle(Particle p, Particle q, bool init){
    double r = distance_rw(p, q);
    double V_LJ = (r < mySys.sigma_cut) ? 4.*mySys.eps*(pow(mySys.sigma/r, 12.)-pow(mySys.sigma/r, 6.)) : 0.;
    return V_LJ;
}


double LJ_potential(bool init){
    double V_LJ_tot = 0.;
    for(int i = 0; i < mySys.NPart; i++){
        for(int j = 0; j < i; j++){
            V_LJ_tot += LJ_potential_per_particle(parts[i], parts[j], init);
        }
        // V_LJ_tot += LJ_potential_per_particle_correction()* mySys.NPart;
    }
    return V_LJ_tot;
}

double kinetic_energy(){
    double en = 0.;
    for(int i = 0; i < mySys.NPart; i++){
        en += 0.5* (parts[i].vx*parts[i].vx + parts[i].vy*parts[i].vy + parts[i].vz*parts[i].vz);
    }
    return en;
}

/*
double compute_pressure(){
    double pressure = 0.;
    double eps = 1.;
    double k_B = 1.;
    double r = 0.;
    double first_virial = 0.;
    double rho = 1.*mySys.NPart/(mySys.box_x*mySys.box_y*mySys.box_z);
    double sigma_cut = mySys.sigma_cut;
    
    for(int i = 0; i < mySys.NPart; i++){
        for(int j = 0; j < i; j++){
            r = distance_rw(parts[i], parts[j]);
            // -r_ij * dV/dr_ij = r_ij * f(r_ij)
            first_virial += 48.*eps*(0.5*pow(mySys.sigma/r, 6.) - pow(mySys.sigma/r, 12)); 
        }
    }
    
    pressure = (mySys.NPart*mySys.T*k_B - first_virial/3.)/(mySys.box_x*mySys.box_y*mySys.box_z);
    
    // pressure tail correction
    pressure += 16./3.*M_PI * (rho*rho) * (2./3. * pow(mySys.sigma/sigma_cut, 9.) - pow(mySys.sigma/sigma_cut, 3.)); 
    
    return pressure;
}
*/

/********************************************************************************************/

double compute_energy_translation(){

    return 0;
}


double compute_energy(int rand_particle, bool init){

    double en = 0.;
    for(int k = 0; k < mySys.NPart; k++){
        if (k == rand_particle) continue;
        en += LJ_potential_per_particle(parts[rand_particle], parts[k], init);
    }
    // en += LJ_potential_per_particle_correction();
   
    return en;

}


double compute_init_energy(){
    return LJ_potential(true);// + kinetic_energy();
}

double compute_init_energy_Harmonic(){
    // Work in 2D here
    if (mySys.exercise == 2 && (mySys.spring_2_x == 0 && mySys.spring_2_y == 0)){ 
        double dx_1 = parts[0].x;
        dx_1 = dx_1 - round(dx_1/mySys.box_x)*mySys.box_x;
        double dy_1 = parts[0].y;
        dy_1 = dy_1 - round(dy_1/mySys.box_y)*mySys.box_y;
        
        return 0.5*mySys.K*(dx_1*dx_1 + dy_1*dy_1); 
    }
    else if (mySys.exercise == 2 && (mySys.spring_2_x != 0 || mySys.spring_2_y != 0)){
        double dx_1 = parts[0].x;
        dx_1 = dx_1 - round(dx_1/mySys.box_x)*mySys.box_x;
        double dy_1 = parts[0].y;
        dy_1 = dy_1 - round(dy_1/mySys.box_y)*mySys.box_y;
        
        double dx_2 = (parts[0].x-mySys.spring_2_x);
        dx_2 = dx_2 - round(dx_2/mySys.box_x)*mySys.box_x;
        double dy_2 = (parts[0].y-mySys.spring_2_y);
        dy_2 = dy_2 - round(dy_2/mySys.box_y)*mySys.box_y;
        return 0.5*mySys.K*(dx_1*dx_1 + dy_1*dy_1) + 0.5*mySys.K*(dx_2*dx_2 + dy_2*dy_2);
    }
    else return 0;
}


/*
Hard Spheres
*/

/*
double compute_init_energy_HS(){
    double res = 0.;
    double dist = 0.;
    for(int i = 0; i < mySys.NPart-1; i++){
        for(int j = i+1; j < mySys.NPart; j++){
            dist = distance_rw(parts[i], parts[j]);
            res += (dist > mySys.sigma) ? 0. : pow(10, 4.);
        }
    }
    return res;
}

double compute_energy_HS(int rand_particle, bool init){
    double en = 0.;
    double r = 0.;
    for(int i = 0; i < mySys.NPart; i++){
        if (rand_particle != i){
            r = (init == true) ? distance_rw(parts[rand_particle], parts[i]) : distance_tmp(parts[rand_particle], parts[i]);
            en += (r > mySys.sigma) ? 0. : pow(10, 4.);
        }
    }
   
    return en;

}
*/

double compute_energy_somemodel(char model[], char type[], int rand_particle, bool init){

    double en = 0.;
    if (model[0] == '0') {
        // Nothing
        en = 0.;    
    } else if (model[0] == '1') {
        // Harmonic potential
        en = compute_init_energy_Harmonic();
    } else {
        printf("Select an energy model!\n");
        exit(1);
    }
    return en;
}
