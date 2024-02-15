// save_data.cpp

#include "save_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>

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


vector<double> readFromFile(const std::string& filePath) {
    vector<double> data;

    // Open the file
    ifstream inputFile(filePath);

    // Check if the file is opened successfully
    if (!inputFile.is_open()) {
        cerr << "Error opening file: " << filePath << endl;
        return data; // Return an empty vector on error
    }

    // Read the file line by line
    string line;
    while (getline(inputFile, line)) {
        // Convert the string to double
        istringstream iss(line);
        double value;
        if (iss >> value) {
            // Add the double value to the vector
            data.push_back(value);
        } else {
            cerr << "Error converting line to double: " << line << endl;
        }
    }

    // Close the file
    inputFile.close();

    return data;
}


void saveVectorToFile(const vector<double>& data, const string& filePath) {
    // Open the file for writing
    ofstream outputFile(filePath);

    // Set precision for output file
    outputFile << setprecision(15);

    // Check if the file is opened successfully
    if (!outputFile.is_open()) {
        cerr << "Error opening file for writing: " << filePath << endl;
        return;
    }

    // Write each element of the vector to the file
    for (const auto& value : data) {
        outputFile << value << endl;
    }

    // Close the file
    outputFile.close();

    cout << "Vector saved to file: " << filePath << endl;
}