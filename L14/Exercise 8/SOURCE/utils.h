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
    
    /*
    // Velocity input
    FILE* fp = fopen(mySys.velocity_file, "r"); 
    if(fp == NULL){ printf("File does not exist!\n"); exit(1); }

    for(int i = 0; i < mySys.NPart; i++) 
    {
        fscanf(fp, "%lf %lf %lf", &a, &b, &c); 
        parts[i].vx = a; parts[i].vy = b; parts[i].vz = c;
    }

    fclose(fp);
    */

    // Check that particles have been initialized correctly
    // printf("%d particles initialized\n", mySys.NPart);
    // for(int i = 0; i < 10; i++){
    //     fprintf(stderr, "x = %lf, y = %lf, z = %lf\n", parts[i].x, parts[i].y, parts[i].z);
    // } 
    // printf("Ok. MC steps...\n");
}


void WriteConf(char filename[], char type[])
{
    FILE* fp = fopen(filename, type); 

    for(int i = 0; i < mySys.NPart; i++) 
    {
        fprintf(fp, "%lf %lf %lf\n", parts[i].x, parts[i].y, parts[i].z); 
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


