void allocate_()
{
    parts = ( Particle *)malloc(mySys.NPart*sizeof(Particle));   
}

void clean_()
{
    free(parts);
}


void ReadConf()
{
    // Position input
    // printf("File Name: '%s'\n", mySys.start_file);
    FILE* fp = fopen(mySys.start_file, "r"); 
    if(fp == NULL){ printf("File does not exist!\n"); exit(1); }

    double a, b, c;

    for(int i = 0; i < mySys.NPart; i++) 
    {
        fscanf(fp, "%lf %lf %lf", &a, &b, &c); 
        parts[i].x = a; parts[i].y = b; parts[i].z = c;
    }

    fclose(fp);
    
    
    // Velocity input
    // printf("File Name: '%s'\n", mySys.velocity_file);
    FILE* fv = fopen("velocity.dat", "r"); 
    if(fv == NULL){ printf("File does not exist!\n"); exit(1); }

    double d, e, f;

    for(int i = 0; i < mySys.NPart; i++) 
    {
        fscanf(fv, "%lf %lf %lf", &d, &e, &f); 
        parts[i].vx = d; parts[i].vy = e; parts[i].vz = f;
    }

    fclose(fv);

    // Check that particles have been initialized correctly
    // printf("%d particles initialized\n", mySys.NPart);
    // for(int i = 0; i < 10; i++){
    //     fprintf(stderr, "x = %lf, y = %lf, z = %lf\n", parts[i].x, parts[i].y, parts[i].z);
    // } 
    // printf("Ok. MC steps...\n");
}


void WriteConf_pos(char filename[], char type[])
{
    FILE* fp = fopen(filename, type); 

    for(int i = 0; i < mySys.NPart; i++) 
    {
        fprintf(fp, "%lf %lf %lf\n", parts[i].x, parts[i].y, parts[i].z); 
    }
    fprintf(fp, "\n");
    fflush(fp); fclose(fp); // Close file when we're done
}

void WriteConf_vel(char filename[], char type[])
{
    FILE* fp = fopen(filename, type); 

    for(int i = 0; i < mySys.NPart; i++) 
    {
        fprintf(fp, "%lf %lf %lf\n", parts[i].vx, parts[i].vy, parts[i].vz); 
    }
    fprintf(fp, "\n");
    fflush(fp); fclose(fp); // Close file when we're done
}



double MinD(double dx, double L){

    double dx1;
    dx1 = dx - rint(dx/L)*L;
    return dx1;
}

double P_Img (double z, double L){
        
    double z1;
    z1 = z - floor(z/L)*L;
    return z1;
}


void BoxMuller(double *x1, double *x2)
{
    double r = sqrt(-2*log(1-drand48()));   // r sampled from p(r)
    double theta = 2.*M_PI*drand48();       // θ sampled from p(θ)
    *x1 = r * cos(theta);
    *x2 = r * sin(theta);
}