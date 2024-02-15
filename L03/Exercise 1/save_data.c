// save_data.c

#include "save_data.h" // Include the corresponding header file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void saveResultsToFile(double results[], int size, char s[100]) {
    FILE *file = fopen(s, "a");
    if (file != NULL) {
        // If the file exists, close it and remove it
        fclose(file);
        remove(s);
    }
    if (file != NULL) {
        file = fopen(s, "a");
        for (int i = 0; i < size; i++) {
            fprintf(file, "%lf\n", results[i]);
        }
        fclose(file);
    }
    else {
        printf("Error opening file for writing!\n");
    }
}

void saveMatrixToFile(double **matrix, int rows, int columns, char s[100]) {
    FILE *file = fopen(s, "w");
    if (file != NULL) {
        // If the file exists, close it and remove it
        fclose(file);
        remove(s);
    }

    if (file != NULL) {
        file = fopen(s, "a");
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                fprintf(file, "%lf ", matrix[i][j]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    } else {
        printf("Error opening file for writing!\n");
    }
}
