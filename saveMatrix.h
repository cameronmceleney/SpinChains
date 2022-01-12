
#ifndef SPINCHAINS_SAVEMATRIX_H
#define SPINCHAINS_SAVEMATRIX_H

#include <fstream>
#include <string>
#include <Eigen3/Eigenvalues>

void saveData(std::string filePath, std::string fileName, Eigen::MatrixXd matrix_to_output)
{
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(filePath+"\\"+fileName);
    if (file.is_open())
    {
        file << matrix_to_output.format(CSVFormat);
        file.close();
    }
}

#endif //SPINCHAINS_SAVEMATRIX_H
