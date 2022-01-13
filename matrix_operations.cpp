#include "matrix_operations.h"

void matrix_operations::saveData(std::string filePath, std::string fileName, Eigen::MatrixXd matrix)
{
    filePath_in = filePath;
    fileName_in = fileName;
    matrix_in = matrix;

    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(filePath_in+"/"+fileName_in);
    if (file.is_open())
    {
        file << matrix_in.format(CSVFormat);
        file.close();
    }
}
