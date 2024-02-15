#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
// #include <cmath>
#include <math.h>
// #include <cmath>

/* our custom libraries */
#include "save_data.h" // Include the custom header file
// #include "RngStream.h"

void display_2(double matrix[2][2]){
    for (int i = 0; i < 2; i++){
        fprintf(stdout, "(");
        for (int j = 0; j < 2; j++){
            fprintf(stdout, "%g", matrix[i][j]);
            if (j != 1){
                fprintf(stdout, ",\t");
            }
            else if (j == 1){
                fprintf(stdout, ")\n");
            }
        }
    }
}

void display_3(double matrix[3][3]){
    for (int i = 0; i < 3; i++){
        fprintf(stdout, "(");
        for (int j = 0; j < 3; j++){
            fprintf(stdout, "%g", matrix[i][j]);
            if (j != 2){
                fprintf(stdout, ",\t");
            }
            else if (j == 2){
                fprintf(stdout, ")\n");
            }
        }
    }
}

int main(){
    /*
    =================
    = Exercise A.5.3 =
    =================
    */
    int rows = 2;
    int cols = 2;
    double P_1[2][2];
    double res[2][2];
    P_1[0][0] = 0.5; P_1[0][1] = 0.5; 
    P_1[1][0] = 1; P_1[1][1] = 0.;
    fprintf(stdout, "P_1 = \n");
    display_2(P_1);

    // Initialize result matrix
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            res[i][j] = P_1[i][j];
        }
    } 

    int N = 100;
    double sum = 0;
    double tmp[2][2];
    for (int n = 0; n < N-1; n++){
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                for (int k = 0; k < 2; k++){
                    sum += P_1[i][k] * res[k][j];
                }
                tmp[i][j] = sum;
                sum = 0;
            }
        }    
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                res[i][j] = tmp[i][j];
            }
        }    
    }
    fprintf(stdout, "\nP_1^%d = \n", N);
    display_2(res);

    // Second part of the exercise:
    double P_2[3][3];
    double result[3][3];
    P_2[0][0] = 0.; P_2[0][1] = 0.; P_2[0][2] = 1.; 
    P_2[1][0] = 0.; P_2[1][1] = 1.; P_2[1][2] = 0.;
    P_2[2][0] = 1./4; P_2[2][1] = 0.; P_2[2][2] = 3./4;
    fprintf(stdout, "\nP_2 = \n");
    display_3(P_2);

    // Initialize result matrix
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            result[i][j] = P_2[i][j];
        }
    } 

    sum = 0;
    double temp[3][3];
    for (int n = 0; n < N-1; n++){
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                for (int k = 0; k < 3; k++){
                    sum += P_2[i][k] * result[k][j];
                }
                temp[i][j] = sum;
                sum = 0;
            }
        }    
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                result[i][j] = temp[i][j];
            }
        }    
    }
    fprintf(stdout, "\nP_2^%d = \n", N);
    display_3(result);


    return 0;
}