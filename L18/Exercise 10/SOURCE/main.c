#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "random.h"
#include "definitions.h"
#include "init.h"
#include "utils.h"
#include "energy.h"
#include "generate_points.h"
#include "VelocityVerlet.h"
#include "MC.h"

int main(int argc, char* argv[] ) 
//start with command line options. char*argv[] means this is an array of indetermined size of pointers to char, in other words an array of strings.argv is a pointer to string array holding the actual parameters.
{

    read_input_file();
 
    allocate_();    			 
    //    or
    // Initialize();  
    do_MC();
    
    // Release memory used by particles, it is not needed anymore
    clean_();

    return 0;
}
// gcc -o ../main main.c -lm

