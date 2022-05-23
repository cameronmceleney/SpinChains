#include <chrono>
#include <ctime>
#include "linspace.h"
#include "SpinChainEigenSolverClass.h"

void SpinChainEigenSolverClass::CalculateEigFreqs() {
    // TODO: rename variables
    // rename all variables following https://manual.gromacs.org/documentation/5.1-current/dev-manual/naming.html

    _totalEquations = GV.GetNumSpins() * 2;

    _fileNameEigenSolver += GV.GetCurrentTime();
    // Compares user's exchange integrals and appends an explanatory string to the filename
    /*if (GV.GetExchangeMinVal() == GV.GetExchangeMaxVal()) {
        _fileNameEigenSolver += "-lin"; // Equal J values means the system has linear exchange ('lin')

    } else if (GV.GetExchangeMinVal() != GV.GetExchangeMaxVal()) {
        _fileNameEigenSolver += "-nonlin"; // Unequal J values means the system has non-linear exchange ('nonlin')

    } else {
        // no statement
    }
    */
    std::cout << "Filename is: " << _fileNameEigenSolver << std::endl; // Informs user of filename to enable directory searching via file explorer search function
    std::cout << "\nFiles will be outputted to: " << GV.GetFilePath() << std::endl; // Showing the selected path will make it easier to find the file

    auto startTimeFindEigens = std::chrono::system_clock::now(); // Separate start time (from system) for the computation of the eigenvalues
    std::time_t startTimeFindEigens_cTimeUse = std::chrono::system_clock::to_time_t( startTimeFindEigens);
    std::cout << "\n------------------------------" << "\nEigenvalues and Eigenvectors";
    std::cout << "\nBegan computation at: " << std::ctime(&startTimeFindEigens_cTimeUse) << std::endl; // Useful for long computations where the start time may be forgotten

    _matrixValues(_totalEquations, _totalEquations); // Generates the matrix but does not allocate memory. That is done as each element is calculated
    _matrixValues = populate_matrix();
    // std::cout << "The populated matrix is :\n" << _matrixValues << std::endl;

    Eigen::EigenSolver <Matrix_xd> eigenSolverMatrixValues(_matrixValues);
    // std::cout << "The eigenvalues are :\n" << eigenSolverMatrixValues.eigenvalues() << std::endl;

    auto stopTimeFindEigens = std::chrono::system_clock::now();
    auto durationTimeFindEigens = std::chrono::duration_cast<std::chrono::milliseconds>(stopTimeFindEigens - startTimeFindEigens);
    std::cout << "Duration to find eigenvectors and values: " << durationTimeFindEigens.count() << "[ms]." << std::endl;

    _matrixValues.resize(0, 0); // Removes large matrix from memory; leads to faster write times for very large matrices

    auto startTimeSaveData = std::chrono::system_clock::now();

    save_data( "eigenvectors_" + _fileNameEigenSolver + ".csv", eigenSolverMatrixValues.eigenvectors().real());
    save_data( "eigenvalues_" + _fileNameEigenSolver + ".csv", eigenSolverMatrixValues.eigenvalues().imag());

    auto stopTimeSaveData = std::chrono::system_clock::now();
    auto durationTimeSaveData = std::chrono::duration_cast<std::chrono::milliseconds>(
            stopTimeSaveData - startTimeSaveData);
    std::cout << "Time to write to files: " << durationTimeSaveData.count() << "[ms]." << std::endl;

    std::time_t stopTimeSaveData_cTimeUse = std::chrono::system_clock::to_time_t(stopTimeSaveData);
    std::cout << "\nFinished computation at: " << std::ctime(&stopTimeSaveData_cTimeUse);
    std::cout << "\n------------------------------\n";

    //std::cout << "Computing V * D * V^(-1) gives: " << std::endl << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

}

void SpinChainEigenSolverClass::PrintVector(std::vector<double> inputVector)
{
    int count = 0;

    std::cout << "size: " << inputVector.size() << std::endl;
    for (double val : inputVector)
    {
        if (++count % 10 == 0)
        {
            std::cout << std::setw(8) << val << " \n";
        }
        else
        {
            std::cout << std::setw(8) << val << " ";
        }
    }
    std::cout << std::endl;
}

void SpinChainEigenSolverClass::save_data( std::string fileName, Matrix_xd generatedMatrix )
{
    /* This class will take a given filename and filepath (supplied by the user), and use that information to create a
     * file for a matrix to be saved into. The matrix must come from the Eigen package and be consistently used throughout
     * all files
     * .
     * This means that if the input (in the argument of the class' object wherever it is called) is of a particular matrix type,
     * then the type of generatedMatrix (line 4 in this file) must be the same as it. An easy way to ensure this is to use
     * Matrix_xd in all instances and then edit the properties of this typedef in the file 'GlobalVariables.h'.*/

    _generatedMatrix = std::move(generatedMatrix);

    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(GV.GetFilePath()+"/"+fileName);

    if (file.is_open()) // Other functions and classes may write/read to this file, so it must be checked that it is opened at this instance
    {
        file << _generatedMatrix.format(CSVFormat); // See documentation for other formats supported by .format
        file.close();
    }
}

Matrix_xd SpinChainEigenSolverClass::populate_matrix()
{

    LinspaceClass exchangeValues{};
    std::vector<double> linspaceExchangeValues; // Holds exchange values for all spins that interact with two other spins

    exchangeValues.set_values(GV.GetExchangeMinVal(), GV.GetExchangeMaxVal(), GV.GetNumSpins()-1, true, false);
    exchangeValues.generate_array();
    // _chainJValues = exchangeValues.build_spinchain();

    /* To simplify the solving of the matrix, setting all unknown frequency variables to zero and then solving the matrix to find eigenvalues proved faster
     * than using an eigen-solver library to find the roots of a characteristic equation populated by angular_frequency values. The outputted
     * eigenvalues of this matrix (matrixToFill) are eigen-frequencies, however they are angular (w). To obtain  frequencies (f), you must find
     * angular_frequency/(2*Pi) for each given eigenvalue. */

    const long totalEquations = GV.GetNumSpins() * 2; // Each spin has an x- and y-axis dependent coupled equation normalised to the z-component. Thus, two equations total

    Matrix_xd matrixToFill(totalEquations, totalEquations); // Each (spin) site location comprises two consecutive rows in matrixToFill

    matrixToFill.setZero(); // Large matrix of known size so more computationally efficient to predefine size in memory

    for (int row = 0, JVal = 0; row < totalEquations; row++) {
        // JVal ensures a spin has the same exchange integral values on its LHS and RHS used in the spin's associated coupled equations
        /* Where an element is matrixToFill(index, index) = 0), this indicates that the element would be a diagonal
         * element of the matrix. While the full computation would be matrixToFill(index, index) =  -1.0  * angular_frequency / _gyroscopicMagneticConstant;),
         * this is an unnecessary series of computations as angular_frequency = 0 is strictly true in this code.*/

        if (row % 2 == 0) {
            // The dm_x/dt coupled equations are the even-numbered rows of the matrix (see notes for details)

            if (row == 0) {
                // Exception for the first dm_x/dt row (1st matrix row) as there is no spin on the LHS of this position and thus no exchange contribution from the LHS
                matrixToFill(row,0) = 0;
                matrixToFill(row,1) = _chainJValues[JVal] + _chainJValues[JVal + 1] + GV.GetBiasField(); //
                matrixToFill(row,3) = -1.0 *  _chainJValues[JVal + 1];
            }
            else if (row > 0 and row < totalEquations - 2) {
                // Handles all other even-numbered rows
                matrixToFill(row,row - 1) = -1.0 *  _chainJValues[JVal];
                matrixToFill(row,row + 0) = 0;
                matrixToFill(row,row + 1) = _chainJValues[JVal] +  _chainJValues[JVal + 1] + GV.GetBiasField();
                matrixToFill(row,row + 3) = -1.0 *  _chainJValues[JVal + 1];
            }
            else if (row == totalEquations - 2) {
                // Exception for the final dm_x/dt row (penultimate matrix row) as there is no spin on the RHS of this position and thus no exchange contribution
                matrixToFill(row,totalEquations - 1) =  _chainJValues[JVal] + _chainJValues[JVal + 1] + GV.GetBiasField(); //
                matrixToFill(row,totalEquations - 2) = 0;
                matrixToFill(row,totalEquations - 3) = -1.0 *  _chainJValues[JVal];
            }
            else {
                // TODO Legacy error handling which needs updating (dm_x/dt rows)
                std::cout << "Error with generating the dx/dt terms on row #{row}. Exiting..." << std::endl;
                std::exit(3);
            }
            continue;
        }
        if (row % 2 == 1) {
            // The dm_y/dt coupled equations are the odd-numbered rows of the matrix (see notes for details)

            if (row == 1) {
                // Exception for the first dm_y/dt row (2nd matrix row) as there is no spin on the LHS of this position and thus no exchange contribution from the LHS
                matrixToFill(row,0) = -1.0 * (_chainJValues[JVal] + _chainJValues[JVal + 1] + GV.GetBiasField()); //
                matrixToFill(row,1) = 0;
                matrixToFill(row,2) =  _chainJValues[JVal + 1];
            }
            else if (row > 1 and row < totalEquations - 1) {
                // Handles all other odd-numbered rows
                matrixToFill(row,row - 3) =  _chainJValues[JVal];
                matrixToFill(row,row - 1) = -1.0 * ( _chainJValues[JVal] +  _chainJValues[JVal + 1] + GV.GetBiasField());
                matrixToFill(row,row + 0) = 0;
                matrixToFill(row,row + 1) =  _chainJValues[JVal + 1];
            }
            else if (row == totalEquations - 1) {
                // Exception for the final dm_y/dt row (final matrix row) as there is no spin on the RHS of this position and thus no exchange contribution
                matrixToFill(row,totalEquations - 1) = 0;
                matrixToFill(row,totalEquations - 2) = -1.0 * ( _chainJValues[JVal] + _chainJValues[JVal + 1] + GV.GetBiasField()); //
                matrixToFill(row,totalEquations - 4) =  _chainJValues[JVal];
            }
            else {
                // TODO Legacy error handling which needs updating (dm_y/dt rows)
                std::cout << "Error with generating the dy/dt terms on row #{row}. Exiting..." << std::endl;
                std::exit(3);
            }
            JVal++;
            continue;
        }
    }

    matrixToFill *= _gyroMagConst; // LLG equation has 'gamma' term outwith the cross-product so all matrix elements must be multiplied by this value

    return matrixToFill;
}