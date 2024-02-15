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


void random_init(){
    int rows = mySys.NPart;
    double low = 0.0;
    double high = mySys.box_x;
    
    double data[rows][SIZE];

    for(int i = 0; i < rows; i++){
        for (int j = 0; j < SIZE; ++j) {
            data[i][j] = rand()/(RAND_MAX*1.)*(high - low) + low;
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