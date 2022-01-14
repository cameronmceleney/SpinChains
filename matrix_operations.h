
#ifndef SPINCHAINS_MATRIX_OPERATIONS_H
#define SPINCHAINS_MATRIX_OPERATIONS_H

#include <fstream>
#include <string>
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

public:
    void save_data(std::string filePath, std::string fileName, Eigen::MatrixXd generatedMatrix );
    Matrix_xd populate_matrix(int numberSpins, double exchangeMinimum, double exchangeMaximum, double biasField, double gyroscopicMagneticConstant);
};


#endif //SPINCHAINS_MATRIX_OPERATIONS_H
