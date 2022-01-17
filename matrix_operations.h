
#ifndef SPINCHAINS_MATRIX_OPERATIONS_H
#define SPINCHAINS_MATRIX_OPERATIONS_H

#include <fstream>
#include <string>
#include <vector>
#include "common.h"

class MatrixOperationsClass {
private:
    std::string _filePath;
    std::string _fileName;
    Eigen::MatrixXd _generatedMatrix;

    int _numberSpins;
    double _exchangeMinimum;
    double _exchangeMaximum;
    double _biasField;
    double _gyroscopicMagneticConstant;
    std::vector<double> _inputVector;

public:
    Matrix_xd populate_matrix(int numberSpins, double exchangeMinimum, double exchangeMaximum, double biasField, double gyroscopicMagneticConstant);
    void save_data(std::string filePath, std::string fileName, Matrix_xd generatedMatrix );
    // Writes a vector to the terminal
    void PrintVector(std::vector<double> inputVector);
};


#endif //SPINCHAINS_MATRIX_OPERATIONS_H
