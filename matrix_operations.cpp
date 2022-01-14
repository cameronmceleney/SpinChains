#include "matrix_operations.h"

#include "linspace.h"

void MatrixOperationsClass::save_data(std::string filePath, std::string fileName, Matrix_xd generatedMatrix )
{

    _filePath = std::move(filePath);
    _fileName = std::move(fileName);
    _generatedMatrix = std::move(generatedMatrix);

    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(_filePath+"/"+_fileName);
    if (file.is_open())
    {
        file << _generatedMatrix.format(CSVFormat);
        file.close();
    }
}

Matrix_xd MatrixOperationsClass::populate_matrix(int numberSpins, double exchangeMinimum, double exchangeMaximum, double biasField, double gyroscopicMagneticConstant)
{

    _numberSpins = numberSpins;
    _exchangeMinimum = exchangeMinimum;
    _exchangeMaximum = exchangeMaximum;
    _biasField = biasField;
    _gyroscopicMagneticConstant = gyroscopicMagneticConstant;

    LinspaceClass Linspace{};
    std::vector<double> JI_values_linspace;
    Linspace.set_values(_exchangeMinimum, _exchangeMaximum, _numberSpins, true);
    JI_values_linspace = Linspace.generate_array(); //calls a linspace function similar to np.linspace()

    static double frequency = 0;
    int linkingExchangesTracker = 0;
    const long totalEquations = _numberSpins * 2;

    //adds a zero to the start and end of JI_values. Insert is faster for large values of numbers compared to push_back()
    std::vector<double> J_array{0}; //initialised with a zero to account for the (P-1)th spin
    J_array.insert(J_array.end(), JI_values_linspace.begin(), JI_values_linspace.end());
    J_array.push_back(0); //add a zero to the end to account for the (N+1)th spin

    Matrix_xd matrixToFill(totalEquations, totalEquations);
    matrixToFill.setZero();

    for (long rowNumber = 0; rowNumber < totalEquations; rowNumber++) {

        if (rowNumber % 2 == 0) {
            //These are even numbered rows, containing the dx/dt terms, so the 'Odd Rule' applies (see notes)

            if (rowNumber == 0) {
                //exception for the first dx/dt row (first row overall), as there is no spin on the LHS of this position
                //static_cast<long double>(-1.0)
                matrixToFill(rowNumber,0) = -1.0  * frequency / _gyroscopicMagneticConstant;
                matrixToFill(rowNumber,1) = (J_array[linkingExchangesTracker] + J_array[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,3) = -1.0 * J_array[linkingExchangesTracker + 1];

            } else if (rowNumber > 0 and rowNumber < totalEquations - 2) {
                //handles all other even numbered rows
                matrixToFill(rowNumber,rowNumber - 1) = -1.0 * J_array[linkingExchangesTracker];
                matrixToFill(rowNumber,rowNumber + 0) = -1.0  * frequency / _gyroscopicMagneticConstant;
                matrixToFill(rowNumber,rowNumber + 1) = (J_array[linkingExchangesTracker] + J_array[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,rowNumber + 3) = -1.0 * J_array[linkingExchangesTracker + 1];

            } else if (rowNumber == totalEquations - 2) {
                // exception for the final dx/dt row (penultimate row overall), as there is no spin on the RHS of this position

                matrixToFill(rowNumber,totalEquations - 1) = (J_array[linkingExchangesTracker] + J_array[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,totalEquations - 2) = -1.0  * frequency / _gyroscopicMagneticConstant;
                matrixToFill(rowNumber,totalEquations - 3) = -1.0 * J_array[linkingExchangesTracker];

            } else {

                std::cout << "Error with generating the dx/dt terms on row #{rowNumber}. Exiting..." << std::endl;
                std::exit(1);

            }
        }

        else if (rowNumber % 2 == 1) {
            //These are odd numbered rows, containing the dy/dt terms, so the 'Even Rule' applies (see notes)
            if (rowNumber == 1) {
                //exception for the first dy/dt row (second row overall), as there is no spin on the LHS of this position
                matrixToFill(rowNumber,0) = -1.0 * (J_array[linkingExchangesTracker] + J_array[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,1) = -1.0  * frequency / _gyroscopicMagneticConstant;
                matrixToFill(rowNumber,2) = J_array[linkingExchangesTracker + 1];

            } else if (rowNumber > 1 and rowNumber < totalEquations - 1) {
                //handles all other odd numbered rows
                matrixToFill(rowNumber,rowNumber - 3) = J_array[linkingExchangesTracker];
                matrixToFill(rowNumber,rowNumber - 1) = -1.0 * (J_array[linkingExchangesTracker] + J_array[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,rowNumber + 0) = -1.0  * frequency / _gyroscopicMagneticConstant;
                matrixToFill(rowNumber,rowNumber + 1) = J_array[linkingExchangesTracker + 1];
            } else if (rowNumber == totalEquations - 1) {
                //exception for the final dy/dt row (final row overall), as there is no spin on the RHS of this position
                matrixToFill(rowNumber,totalEquations - 1) = -1.0  * frequency / _gyroscopicMagneticConstant;
                matrixToFill(rowNumber,totalEquations - 2) = -1.0 * (J_array[linkingExchangesTracker] + J_array[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,totalEquations - 4) = J_array[linkingExchangesTracker];
            } else {

                std::cout << "Error with generating the dy/dt terms on row #{rowNumber}. Exiting..." << std::endl;
                std::exit(1);
            }
            linkingExchangesTracker += 1;
        }

        else {
            std::cout << "An issue arose in generating the matrix. Exiting..." << std::endl;
            std::exit(1);
        }
    }

    matrixToFill *= _gyroscopicMagneticConstant;
    return matrixToFill;
}
