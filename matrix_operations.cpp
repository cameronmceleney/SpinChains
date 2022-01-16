#include "matrix_operations.h"
#include "linspace.h"

void MatrixOperationsClass::PrintVector(std::vector<double> inputVector)
{
    _inputVector = inputVector;
    std::cout << "size: " << _inputVector.size() << std::endl;
    for (double val : _inputVector)
        std::cout << val << " ";
    std::cout << std::endl;
}

void MatrixOperationsClass::save_data(std::string filePath, std::string fileName, Matrix_xd generatedMatrix )
{
    /* This class will take a given filename and filepath (supplied by the user), and use that information to create a
     * file for a matrix to be saved into. The matrix must come from the Eigen package and be consistently used throughout
     * all files
     * .
     * This means that if the input (in the argument of the class' object wherever it is called) is of a particular matrix type,
     * then the type of generatedMatrix (line 4 in this file) must be the same as it. An easy way to ensure this is to use
     * Matrix_xd in all instances and then edit the properties of this typedef in the file 'common.h'.*/

    _filePath = std::move(filePath);
    _fileName = std::move(fileName);
    _generatedMatrix = std::move(generatedMatrix);

    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(_filePath+"/"+_fileName);

    if (file.is_open()) // Other functions and classes may write/read to this file, so it must be checked that it is opened at this instance
    {
        file << _generatedMatrix.format(CSVFormat); // See documentation for other formats supported by .format
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
    std::vector<double> linspaceExchangeValues; // Holds exchange values for all spins that interact with two other spins

    Linspace.set_values(_exchangeMinimum, _exchangeMaximum, _numberSpins, true);
    linspaceExchangeValues = Linspace.generate_array();

    /* To simplify the solving of the matrix, setting all unknown frequency variables to zero and then solving the matrix to find eigenvalues proved faster
     * than using an eigen-solver library to find the roots of a characteristic equation populated by angular_frequency values. The outputted
     * eigenvalues of this matrix (matrixToFill) are eigen-frequencies, however they are angular (w). To obtain  frequencies (f), you must find
     * angular_frequency/(2*Pi) for each given eigenvalue. */

    // static double angular_frequency = 0;
    int linkingExchangesTracker = 0; // Ensures a spin has the same exchange integral values on its LHS and RHS used in the spin's associated coupled equations
    const long totalEquations = _numberSpins * 2; // Each spin has an x- and y-axis dependent coupled equation normalised to the z-component. Thus, two equations total

    // Insert is faster for large values of numbers compared to push_back().
    std::vector<double> fullChainExchangeValues{0}; // Initialised with a zero to account for the exchange from the (P-1)th LHS spin
    fullChainExchangeValues.insert(fullChainExchangeValues.end(), linspaceExchangeValues.begin(), linspaceExchangeValues.end());
    fullChainExchangeValues.push_back(0); // Appends a zero to the end to account for the exchange from the (N+1)th RHS spin

    Matrix_xd matrixToFill(totalEquations, totalEquations); // Each (spin) site location comprises two consecutive rows in matrixToFill

    matrixToFill.setZero(); // Large matrix of known size so more computationally efficient to predefine size in memory

    for (long rowNumber = 0; rowNumber < totalEquations; rowNumber++) {
        /* Where an element is matrixToFill(index, index) = 0), this indicates that the element would be a diagonal
         * element of the matrix. While the full computation would be matrixToFill(index, index) =  -1.0  * angular_frequency / _gyroscopicMagneticConstant;),
         * this is an unnecessary series of computations as angular_frequency = 0 is strictly true in this code.*/

        if (rowNumber % 2 == 0) {
            // The dm_x/dt coupled equations are the even-numbered rows of the matrix (see notes for details)

            if (rowNumber == 0) {
                // Exception for the first dm_x/dt row (1st matrix row) as there is no spin on the LHS of this position and thus no exchange contribution from the LHS
                matrixToFill(rowNumber,0) = 0;
                matrixToFill(rowNumber,1) = (fullChainExchangeValues[linkingExchangesTracker] + fullChainExchangeValues[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,3) = -1.0 * fullChainExchangeValues[linkingExchangesTracker + 1];

            } else if (rowNumber > 0 and rowNumber < totalEquations - 2) {
                // Handles all other even-numbered rows
                matrixToFill(rowNumber,rowNumber - 1) = -1.0 * fullChainExchangeValues[linkingExchangesTracker];
                matrixToFill(rowNumber,rowNumber + 0) = 0;
                matrixToFill(rowNumber,rowNumber + 1) = (fullChainExchangeValues[linkingExchangesTracker] + fullChainExchangeValues[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,rowNumber + 3) = -1.0 * fullChainExchangeValues[linkingExchangesTracker + 1];

            } else if (rowNumber == totalEquations - 2) {
                // Exception for the final dm_x/dt row (penultimate matrix row) as there is no spin on the RHS of this position and thus no exchange contribution
                matrixToFill(rowNumber,totalEquations - 1) = (fullChainExchangeValues[linkingExchangesTracker] + fullChainExchangeValues[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,totalEquations - 2) = 0;
                matrixToFill(rowNumber,totalEquations - 3) = -1.0 * fullChainExchangeValues[linkingExchangesTracker];

            } else {
                // TODO Legacy error handling which needs updating (dm_x/dt rows)
                std::cout << "Error with generating the dx/dt terms on row #{rowNumber}. Exiting..." << std::endl;
                std::exit(1);

            }
        }

        else if (rowNumber % 2 == 1) {
            // The dm_y/dt coupled equations are the odd-numbered rows of the matrix (see notes for details)

            if (rowNumber == 1) {
                // Exception for the first dm_y/dt row (2nd matrix row) as there is no spin on the LHS of this position and thus no exchange contribution from the LHS
                matrixToFill(rowNumber,0) = -1.0 * (fullChainExchangeValues[linkingExchangesTracker] + fullChainExchangeValues[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,1) = 0;
                matrixToFill(rowNumber,2) = fullChainExchangeValues[linkingExchangesTracker + 1];

            } else if (rowNumber > 1 and rowNumber < totalEquations - 1) {
                // Handles all other odd-numbered rows
                matrixToFill(rowNumber,rowNumber - 3) = fullChainExchangeValues[linkingExchangesTracker];
                matrixToFill(rowNumber,rowNumber - 1) = -1.0 * (fullChainExchangeValues[linkingExchangesTracker] + fullChainExchangeValues[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,rowNumber + 0) = 0;
                matrixToFill(rowNumber,rowNumber + 1) = fullChainExchangeValues[linkingExchangesTracker + 1];

            } else if (rowNumber == totalEquations - 1) {
                // Exception for the final dm_y/dt row (final matrix row) as there is no spin on the RHS of this position and thus no exchange contribution
                matrixToFill(rowNumber,totalEquations - 1) = 0;
                matrixToFill(rowNumber,totalEquations - 2) = -1.0 * (fullChainExchangeValues[linkingExchangesTracker] + fullChainExchangeValues[linkingExchangesTracker + 1] + _biasField);
                matrixToFill(rowNumber,totalEquations - 4) = fullChainExchangeValues[linkingExchangesTracker];

            } else {
                // TODO Legacy error handling which needs updating (dm_y/dt rows)
                std::cout << "Error with generating the dy/dt terms on row #{rowNumber}. Exiting..." << std::endl;
                std::exit(1);
            }

            linkingExchangesTracker += 1; // Solving one set of coupled equations means the tracker can be increased to reflect moving along one spin site position in the chain
        }

        else {
            // TODO Legacy error handling which needs updating (matrixToFill)
            std::cout << "An issue arose in generating the matrix. Exiting..." << std::endl;
            std::exit(1);
        }
    }

    matrixToFill *= _gyroscopicMagneticConstant; // LLG equation has 'gamma' term outwith the cross-product so all matrix elements must be multiplied by this value

    return matrixToFill;
}
