
#ifndef SPINCHAINS_MATRIX_OPERATIONS_H
#define SPINCHAINS_MATRIX_OPERATIONS_H

#include <fstream>
#include <string>
#include "common.h"

class matrix_operations {
private:
    std::string filePath_in;
    std::string fileName_in;
    Eigen::MatrixXd matrix_in;

    int NumSpins_in;
    double ExchangeIntegral_min_in;
    double ExchangeIntegral_max_in;
    double BiasField_in;
    double GyroMagConst_in;

public:
    void SaveData(std::string filePath, std::string fileName, Eigen::MatrixXd matrix_to_output );
    MatrixcXd PopulateMatrix(int NumSpins, double ExchangeIntegral_min, double ExchangeIntegral_max, double BiasField, double GyroMagConst);
};


#endif //SPINCHAINS_MATRIX_OPERATIONS_H
