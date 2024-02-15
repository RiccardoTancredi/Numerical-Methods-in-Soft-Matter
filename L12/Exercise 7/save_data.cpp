// save_data.cpp

#include "save_data.h"
#include <iostream>
#include <fstream>

using namespace std;

void saveResultsToFile(double results[], int size, const string& filename) {
    ofstream file(filename, ios::app); // Open the file in append mode

    if (file.is_open()) {
        // If the file exists, close it and remove it
        file.close();
        remove(filename.c_str());
    }

    file.open(filename, ios::app); // Open the file again in append mode

    if (file.is_open()) {
        for (int i = 0; i < size; i++) {
            file << results[i] << endl;
        }
        file.close();
    } else {
        cout << "Error opening file for writing!" << endl;
    }
}

void saveMatrixToFile(double** matrix, int rows, int columns, const string& filename) {
    ofstream file(filename, ios::out); // Open the file in write mode

    if (file.is_open()) {
        // If the file exists, close it and remove it
        file.close();
        remove(filename.c_str());
    }

    file.open(filename, ios::app); // Open the file again in append mode

    if (file.is_open()) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                file << matrix[i][j] << " ";
            }
            file << endl;
        }
        file.close();
    } else {
        cout << "Error opening file for writing!" << endl;
    }
}

void saveMatrixToFile_int(int** matrix, int rows, int columns, const string& filename) {
    ofstream file(filename, ios::out); // Open the file in write mode

    if (file.is_open()) {
        // If the file exists, close it and remove it
        file.close();
        remove(filename.c_str());
    }

    file.open(filename, ios::app); // Open the file again in append mode

    if (file.is_open()) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                file << matrix[i][j] << " ";
            }
            file << endl;
        }
        file.close();
    } else {
        cout << "Error opening file for writing!" << endl;
    }
}
