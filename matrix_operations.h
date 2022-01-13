
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

public:
    void saveData(std::string filePath, std::string fileName, Eigen::MatrixXd matrix_to_output );
};


#endif //SPINCHAINS_MATRIX_OPERATIONS_H
