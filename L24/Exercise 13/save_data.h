// save_data.h

#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <string>
#include <vector>

using namespace std;

void saveResultsToFile(double results[], int size, const string& filename);
void saveMatrixToFile(double** matrix, int rows, int columns, const string& filename);
void saveMatrixToFile_int(int** matrix, int rows, int columns, const string& filename);
vector<double> readFromFile(const string& filePath);
void saveVectorToFile(const vector<double>& data, const string& filePath);

#endif
