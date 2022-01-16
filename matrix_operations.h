
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
    Matrix_xd populate_matrix(int numberSpins, double exchangeMinimum, double exchangeMaximum, double biasField, double gyroscopicMagneticConstant);
    void save_data(std::string filePath, std::string fileName, Matrix_xd generatedMatrix );
};


#endif //SPINCHAINS_MATRIX_OPERATIONS_H
