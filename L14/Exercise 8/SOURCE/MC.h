void copycoordinates(int idp)
{
    parts[idp].ox = parts[idp].x;
    parts[idp].oy = parts[idp].y;
    parts[idp].oz = parts[idp].z;
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

void do_MC(){

    srand48(mySys.seed);

    char folder_name[1000] = "";
    char file_path_1[1000] = "";
    char file_path_2[1000] = "";
    char file_path_3[1000] = "";
    char restartname[1000] = "";

    char temp_kind[10]; 
    if (mySys.kind[0] == 'r'){
        strcpy(temp_kind, "random");
    }
    else{
        strcpy(temp_kind, "cube");
    }

    if (mySys.T == 1) {
        if (fabs(mySys.disp_max - 0.01) < 1e-6 || fabs(mySys.disp_max - 0.05) < 1e-6) {
            snprintf(folder_name, sizeof(folder_name), "res/%s/%d_%.2f/%d/", temp_kind, (int)mySys.box_x, mySys.disp_max, mySys.iter);
        } else if (fabs(mySys.disp_max - 0.3) < 1e-6 || fabs(mySys.disp_max - 0.5) < 1e-6) {
            snprintf(folder_name, sizeof(folder_name), "res/%s/%d_%.1f/%d/", temp_kind, (int)mySys.box_x, mySys.disp_max, mySys.iter);
        } else {
            snprintf(folder_name, sizeof(folder_name), "res/%s/%d_%d/%d/", temp_kind, (int)mySys.box_x, (int)mySys.disp_max, mySys.iter);
        }
    } else {

        if (mySys.T < 1) {
            snprintf(folder_name, sizeof(folder_name), "res/pressure/T_%.1f/", mySys.T);
        } else {
            snprintf(folder_name, sizeof(folder_name), "res/pressure/T_%d/", (int)mySys.T);
        }
    }

    // Create data:
    if (temp_kind[0] == 'c') {
        generate_init_configuration();
        printf("File created!\n");
    }
    else if (mySys.T != 1){
        random_init();
        printf("File created!\n");
    }
    
    ReadConf();

    snprintf(file_path_1, sizeof(file_path_1), "%sV_%d_energy.dat", folder_name, (int)mySys.box_x);
    snprintf(file_path_2, sizeof(file_path_2), "%sV_%d_acceptance.dat", folder_name, (int)mySys.box_x);
    snprintf(file_path_3, sizeof(file_path_3), "%sV_%d_pressure.dat", folder_name, (int)mySys.box_x);
    snprintf(restartname, sizeof(restartname), "%sV_%d_restartpoint.dat", folder_name, (int)mySys.box_x);

    printf("Folder name: %s\n", file_path_1);

    FILE* f = fopen(file_path_1, "w");
    FILE* g = fopen(file_path_2, "w");
    FILE* p = fopen(file_path_3, "w"); 

    char model[2];
    if(mySys.T == 1) strcpy(model, "1");    // H-S
    else strcpy(model, "0");    // L-J
    
    // double energy = compute_init_energy(); // initial energy
    double energy = compute_energy_somemodel(model, "init", 0, 0);
    printf("\nThe initial energy is E_in = %lf\n", energy);
    mySys.accept = 0;
    for(mySys.step=0; mySys.step < mySys.NSteps; mySys.step++){
        energy = do_MC_sweep(energy, model);
        if(mySys.step > mySys.NPrint) fprintf(f, "%.15f\n", energy); 
        if(mySys.step > mySys.NPrint) fprintf(g, "%d\n", mySys.accept);
        if(model[0] == '0' && mySys.step > mySys.NPrint) fprintf(p, "%.15f\n", compute_pressure());

        // if(mySys.step % mySys.NPrint == 0){
        //     if (mySys.step == 0) WriteConf(restartname, "w");
        //     else WriteConf(restartname, "a");
        // }
     
        // if(mySys.step % mySys.NPrint == 0){ 
        //     printf("dumping...\n");
        // }
    }
    if(model[0] == '0') printf("Final pressure = %lf\n", compute_pressure());
   
    fclose(f); fclose(g);
   
}

