void copycoordinates(int idp)
{
    parts[idp].ox = parts[idp].x;
    parts[idp].oy = parts[idp].y;
    parts[idp].oz = parts[idp].z;
}

void copyvelocities(int idp)
{
    parts[idp].ovx = parts[idp].vx;
    parts[idp].ovy = parts[idp].vy;
    parts[idp].ovz = parts[idp].vz;
}

void rollback(int idp)
{
    parts[idp].x = parts[idp].ox;
    parts[idp].y = parts[idp].oy;
    parts[idp].z = parts[idp].oz;
}


void translation(int idp)
{
    parts[idp].x = parts[idp].tx;
    parts[idp].y = parts[idp].ty;
    parts[idp].z = parts[idp].tz;
}

void Velocity_rescaling(double k_energy){
    double T_k = k_energy*2./(3.*mySys.NPart); 
    double lambda = sqrt(mySys.T/T_k);
    for(int o = 0; o < mySys.NPart; o++){
        parts[o].vx *= lambda; 
        parts[o].vy *= lambda; 
        parts[o].vz *= lambda; 
    }
}

void Andersen_thermostat(double prob){
    double r_i = 0.;
    double r_gauss[3];
    double theta_gauss[3];
    for(int o = 0; o < mySys.NPart; o++){
        r_i = ((double)rand() / RAND_MAX*1.);
        if (r_i < prob){   
            // Sample velocities from a Normal distribution N(0, sqrt(T/m)=1) 
            for(int dim = 0; dim < 3; dim++){
                r_gauss[dim] = sqrt(-2.*log(1.-((double)rand() / RAND_MAX*1.)))*sqrt(mySys.T); // r sampled from p(r), scaled with the variace T/m of the Maxwell-Boltzmann distribution
                theta_gauss[dim] = 2.*M_PI*((double)rand() / RAND_MAX*1.); // θ sampled from p(θ)
            }
            parts[o].vx = r_gauss[0]*cos(theta_gauss[0]); 
            parts[o].vy = r_gauss[1]*cos(theta_gauss[1]); 
            parts[o].vz = r_gauss[2]*cos(theta_gauss[2]); 
        }
    }
}

void PBC(int idp)
{
    //Spherical coordinates
	double rand_r = rand()/(RAND_MAX*1.)*mySys.disp_max; 
	double rand_phi = rand()/(RAND_MAX*1.)*2.*M_PI;
	double rand_theta = rand()/(RAND_MAX*1.)*M_PI;
	
	//Switch to cartesian coordinates
    double rand_num_x = rand_r*sin(rand_theta)*cos(rand_phi);
    double rand_num_y = rand_r*sin(rand_theta)*sin(rand_phi);
    double rand_num_z = rand_r*cos(rand_theta);

    parts[idp].tx = P_Img(parts[idp].x + rand_num_x, mySys.box_x);
    parts[idp].ty = P_Img(parts[idp].y + rand_num_y, mySys.box_y);
    parts[idp].tz = P_Img(parts[idp].z + rand_num_z, mySys.box_z);
}



double do_MC_sweep(double energy, char model[])
{
    for(int k = 0; k < mySys.NPart; k++){

        // 1. Select a particle at random: set seed and generate a random number 
        // double randomValue = (double)rand() / (RAND_MAX*1.); 
        int rand_particle = floor(rand()/(RAND_MAX*1.) * mySys.NPart);
        
        // 2. Propose a displacement in each direction. The maximum displacement should be set as a parameter d_max.
        // 3. Compute the energy of the system before the displacement -> init_energy
        double init_energy = compute_energy_somemodel(model, "other", rand_particle, true);
        // 4. Displace the particle and compute the energy of the system after the displacement    
        // Periodic Boundary Conditions
        PBC(rand_particle);
        double fin_energy = compute_energy_somemodel(model, "other", rand_particle, false);
        double Delta_E = fin_energy - init_energy;
        // 5. Accept or reject according to the Metropolis rule.
        
        // Hard Spheres:
        if (model[0] == '1'){     
            if (Delta_E < pow(10, 4.)){
                // accept the move
                mySys.accept += 1;

                // // Set old coordinates
                // copycoordinates(rand_particle);
                
                // Set new accepted coordinates
                translation(rand_particle);
                
                // Update energy
                energy += Delta_E;
            } 
        }
        
        else if(model[0] == '0'){
            if(Delta_E < 0 || log((double)rand() / RAND_MAX*1.) < -Delta_E/(mySys.T)){
                // accept the move
                mySys.accept += 1;

                // // Set old coordinates
                // copycoordinates(rand_particle);
                
                // Set new accepted coordinates
                translation(rand_particle);

                // Update energy
                energy += Delta_E;
            } 
        }
    }
    return energy;
}

void do_move(double dt){
    double r_ij = 0;
    double dxu = 0, dyu = 0, dzu = 0;  
        
    for(int i = 0; i < mySys.NPart; i++){

        // double x_t = parts[i].ox; double y_t = parts[i].oy; double z_t = parts[i].oz;
        // double p_x_t = parts[i].ovx; double p_y_t = parts[i].ovy; double p_z_t = parts[i].ovz;                
        
        set_forces(i);  // set forces of particel i to 0
        double force = 0.;
        for(int j = 0; j < mySys.NPart; j++){
            if(i == j) continue;
            r_ij = distance_rw(parts[i], parts[j]);
            dxu = parts[i].ox - parts[j].ox;
            dyu = parts[i].oy - parts[j].oy;
            dzu = parts[i].oz - parts[j].oz;

            // Periodic Boundary Conditions
            dxu = dxu - round(dxu/mySys.box_x)*mySys.box_x;
            dyu = dyu - round(dyu/mySys.box_y)*mySys.box_y;
            dzu = dzu - round(dzu/mySys.box_z)*mySys.box_z;

            force = compute_force(r_ij);
            parts[i].fx += force*dxu; 
            parts[i].fy += force*dyu;
            parts[i].fz += force*dzu;
        }
    }

    // Update particles' position and velocity using the
    // Velocity Verlet algortihm:
    Velocity_Verlet(dt);
}

void do_MC(){

    srand48(mySys.seed);

    char folder_name[1000] = "";
    char file_path_1[1000] = "";
    char file_path_2[1000] = "";
    char file_path_3[1000] = "";
    char file_path_4[1000] = "";
    
    // No thermostats
    if (mySys.thermostat == 0){    
        snprintf(folder_name, sizeof(folder_name), "res/Velocity_Verlet/%.1f/", mySys.sigma_cut);
    }
    // Velocity rescaling: non canonical thermostat
    else if(mySys.thermostat == 1){
        snprintf(folder_name, sizeof(folder_name), "res/Thermostats/V_rescaling/%d/", mySys.NPart);
    }
    // Andersen thermostat: canonical thermostat
    else{
        snprintf(folder_name, sizeof(folder_name), "res/Thermostats/Andersen/%d/", mySys.NPart);
    }

    snprintf(file_path_1, sizeof(file_path_1), "%senergy.dat", folder_name);
    snprintf(file_path_2, sizeof(file_path_2), "%sdata.dat", folder_name);
    snprintf(file_path_3, sizeof(file_path_3), "%svelocity.dat", folder_name);
    snprintf(file_path_4, sizeof(file_path_4), "%skinetic_energy.dat", folder_name);
    
    printf("Folder name: %s\n", file_path_1);
    
    // Read data:    
    ReadConf();

    FILE* f = fopen(file_path_1, "w");
    FILE* g = fopen(file_path_4, "w");
    // FILE* h = fopen(file_path_3, "w");

    double dt = pow(10, -3.);

    char model[2] = "0";    // L-J
    // if(mySys.T == 1) strcpy(model, "1");    // H-S
    // else strcpy(model, "0");    
    
    // Copy coordinates
    for(int l = 0; l < mySys.NPart; l++){
        copycoordinates(l);
        copyvelocities(l);
    }

    if(mySys.thermostat == 0) printf("\nsigma_cut = %.1f\n", mySys.sigma_cut);

    // initial energy
    double energy = compute_energy_somemodel(model, "init", 0, 0);
    printf("\nThe initial energy is E_in = %lf\n", energy);
    double k_energy = kinetic_energy();
    mySys.accept = 0;
    double omega = 2.*pow(10, 3.);
    
    for(mySys.step=0; mySys.step < mySys.NSteps; mySys.step++){
        
        do_move(dt);

        // Velocity rescaling
        if(mySys.thermostat == 1){
            Velocity_rescaling(kinetic_energy());
        }
        // Andersen thermostat
        else if (mySys.thermostat == 2)
        {
            Andersen_thermostat(omega*dt);
        }
        
        
        // Copy coordinates
        for(int l = 0; l < mySys.NPart; l++){
            copycoordinates(l);
            copyvelocities(l);
        }

        energy = compute_energy_somemodel(model, "init", 0, 0);
        k_energy = kinetic_energy();
        fprintf(f, "%.15f\n", energy);
        fprintf(g, "%.15f\n", k_energy); 

        // if(mySys.step > mySys.NPrint) fprintf(g, "%d\n", mySys.accept);
        // if(model[0] == '0' && mySys.step > mySys.NPrint) fprintf(p, "%.15f\n", compute_pressure());

        // if(mySys.step % mySys.NPrint == 0){
        if (mySys.step == 0) {
            WriteConf_pos(file_path_2, "w");
            WriteConf_vel(file_path_3, "w");
        }
        else {
            WriteConf_pos(file_path_2, "a");
            WriteConf_vel(file_path_3, "a");
        }
        // }

    }
    // if(model[0] == '0') printf("Final pressure = %lf\n", compute_pressure());
    printf("Final energy = %.15f\n", energy);
    fclose(f); // fclose(g);
   
}

