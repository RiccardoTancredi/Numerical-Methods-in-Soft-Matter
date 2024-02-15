void read_input_file()
{
    FILE *_myfile = fopen("param.dat", "r");
    
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
   
      if(strcmp("N",_key) == 0){ mySys.NPart = atoi(_value); continue; }
      if(strcmp("N_steps",_key) == 0){ mySys.NSteps = atoi(_value); continue; }
      if(strcmp("box_x",_key) == 0){ mySys.box_x = atof(_value); continue; } 
      if(strcmp("box_y",_key) == 0){ mySys.box_y = atof(_value); continue; }
      if(strcmp("box_z",_key) == 0){ mySys.box_z = atof(_value); continue; }
      if(strcmp("temperature",_key) == 0){ mySys.T = atof(_value); continue; }
      if(strcmp("delta_x",_key) == 0){ mySys.disp_max = atof(_value); continue; }
      if(strcmp("myseed",_key) == 0){ mySys.seed = atoi(_value); continue; } 
      if(strcmp("sigma",_key) == 0){ mySys.sigma = atof(_value); continue; }
      //if(strcmp("model",_key) == 0){ mySys.model = atoi(_value); continue; }
      if(strcmp("ratio_samplings",_key) == 0){ mySys.NPrint = atoi(_value); continue; }
      if(strcmp("iteration",_key) == 0){ mySys.iter = atoi(_value); continue; }    
      if (strcmp("start_file", _key) == 0) {
        strncpy(mySys.start_file, _value, sizeof(mySys.start_file) - 1);
        mySys.start_file[sizeof(mySys.start_file) - 1] = '\0';  // Ensure null-terminated string

        // Remove trailing newline if it exists
        size_t length = strlen(mySys.start_file);
        if (length > 0 && (mySys.start_file[length - 1] == '\n' || mySys.start_file[length - 1] == '\r')) {
            mySys.start_file[length - 1] = '\0';
        }
        continue;
      }
      if (strcmp("kind", _key) == 0) {
        strncpy(mySys.kind, _value, sizeof(mySys.kind) - 1);
        mySys.kind[sizeof(mySys.kind) - 1] = '\0';  // Ensure null-terminated string

        // Remove trailing newline if it exists
        size_t length_kind = sizeof(mySys.kind);
        if (length_kind > 0 && (mySys.kind[length_kind - 1] == '\n' || mySys.kind[length_kind - 1] == '\r')) {
            mySys.kind[length_kind - 1] = '\0';
        }
        continue;
      }

      // if (strcmp("velocity_file", _key) == 0) {
      //   strncpy(mySys.velocity_file, _value, sizeof(mySys.velocity_file) - 1);
      //   mySys.velocity_file[sizeof(mySys.velocity_file) - 1] = '\0';  // Ensure null-terminated string

      //   // Remove trailing newline if it exists
      //   size_t length_v = strlen(mySys.velocity_file);
      //   if (length_v > 0 && mySys.velocity_file[length_v - 1] == '\n') {
      //       mySys.velocity_file[length_v - 1] = '\0';
      //   } 

      //   continue;
      // }
      
    }
    
    fclose(_myfile);

}

