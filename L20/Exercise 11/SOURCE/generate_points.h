#define SIZE 3  // Assuming 3D points

void generate_init_configuration() {
    int rows = mySys.NPart;
    double lattice_side = mySys.box_x;
    
    double data[rows][SIZE];

    // Seed the random number generator
    srand48((long)mySys.seed);

    int counter = 0; 
	for(int i = 0; i < lattice_side; i++)
	{
		for(int j = 0; j < lattice_side; j++)
		{
			for(int k = 0; k < lattice_side; k++)
			{
                data[counter][0] = i;
                data[counter][1] = j; 
                data[counter][2] = k;
                counter +=1;
				if(counter >= mySys.NPart) break;
			}
			if(counter >= mySys.NPart) break;
		}
		if(counter >= mySys.NPart) break;
	}
	
	return;
    
    // Save the generated points
    FILE *fp = fopen(mySys.start_file, "w");

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            fprintf(fp, "%lf ", data[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

}

double distance(double *p, double *q){
    // compute distance
    double dxu = p[0] - q[0];
    double dyu = p[1] - q[1];
    double dzu = p[2] - q[2];
    double dx = dxu - round(dxu/mySys.box_x)*mySys.box_x;
    double dy = dyu - round(dyu/mySys.box_y)*mySys.box_y;
    double dz = dzu - round(dzu/mySys.box_z)*mySys.box_z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}


void random_init(){
    int rows = mySys.NPart;
    double low = 0.0;
    double high = mySys.box_x;
    
    double data[rows][SIZE];
    data[0][0] = rand()/(RAND_MAX*1.)*(high - low) + low;
    double trial[1][SIZE];
    bool go;
    int counter = 0;
    for(int i = 1; i < rows; i++){
        go = true;
        while(go){
            counter = 0;
            for (int j = 0; j < SIZE; ++j) {
                trial[0][j] = rand()/(RAND_MAX*1.)*(high - low) + low;
		    }
            // compute distances
            for(int k = 0; k < i; k++){
                if(distance(data[k], trial[0]) > mySys.sigma_cut) counter +=1;
            }    
            if(counter == i) {
                for (int j = 0; j < SIZE; ++j) {
                    data[i][j] = trial[0][j];
                }
                go = false;
            }
        }
	}
    
    // Save the generated points
    FILE *fp = fopen(mySys.start_file, "w");

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            fprintf(fp, "%lf ", data[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

}