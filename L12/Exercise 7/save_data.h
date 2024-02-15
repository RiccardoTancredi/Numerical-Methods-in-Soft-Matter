// save_data.h

#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <string>
using namespace std;

void saveResultsToFile(double results[], int size, const string& filename);
void saveMatrixToFile(double** matrix, int rows, int columns, const string& filename);
void saveMatrixToFile_int(int** matrix, int rows, int columns, const string& filename);

#endif
