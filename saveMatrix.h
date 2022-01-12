
#ifndef SPINCHAINS_SAVEMATRIX_H
#define SPINCHAINS_SAVEMATRIX_H

#include <fstream>
#include <string>
#include <Eigen3/Eigenvalues>

class saveMatrix {
private:
    std::string filePath_in;
    std::string fileName_in;
    Eigen::MatrixXd matrix_in;

public:
    void saveData(std::string filePath, std::string fileName, Eigen::MatrixXd matrix_to_output );
};


#endif //SPINCHAINS_SAVEMATRIX_H
