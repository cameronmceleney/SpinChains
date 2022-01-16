#include <chrono>
#include <cmath>
#include <ctime>
#include "matrix_operations.h"
#include "SpinChainEigenSolverClass.h"

void SpinChainEigenSolverClass::SolveInputs(int numberOfSpins) {

    _numberOfSpins = numberOfSpins;

    // TODO: rename variables
    // rename all variables following https://manual.gromacs.org/documentation/5.1-current/dev-manual/naming.html

    _fileLocation = R"(/Users/cameronmceleney/OneDrive - University of Glasgow/University/PhD/1st Year/C++/Chainspins/Data Outputs)";

    _totalEquations = _numberOfSpins * 2;

    _fileNameNumSpins = std::to_string(_numberOfSpins) + "spins";

    std::cout << "Enter the minimum exchange value: ";
    std::cin >> _exchangeMinimum;

    std::cout << "Enter the maximum exchange value: ";
    std::cin >> _exchangeMaximum;

    // Compares user's exchange integrals and appends an explanatory string to the filename
    if (_exchangeMinimum == _exchangeMaximum) {
        _fileNameNumSpins += "-lin"; // Equal J values means the system has linear exchange ('lin')

    } else if (_exchangeMinimum != _exchangeMaximum) {
        _fileNameNumSpins += "-nonlin"; // Unequal J values means the system has non-linear exchange ('nonlin')

    } else {
        /* no statement */
    }

    MatrixOperationsClass MatrixOperations;

    std::cout << "Filename is: " << _fileNameNumSpins << std::endl; // Informs user of filename to enable directory searching via file explorer search function
    std::cout << "\nFiles will be outputted to: " << _fileLocation << std::endl; // Showing the selected path will make it easier to find the file

    auto startTimeFindEigens = std::chrono::system_clock::now(); // Separate start time (from system) for the computation of the eigenvalues
    std::time_t startTimeFindEigens_cTimeUse = std::chrono::system_clock::to_time_t( startTimeFindEigens);
    std::cout << "\nBegan computation at: " << std::ctime(&startTimeFindEigens_cTimeUse) << std::endl; // Useful for long computations where the start time may be forgotten

    _matrixValues(_totalEquations, _totalEquations); // Generates the matrix but does not allocate memory. That is done as each element is calculated
    _matrixValues = MatrixOperations.populate_matrix(_numberOfSpins, _exchangeMinimum, _exchangeMaximum, 0.1,29.2 * 2 * M_PI);
    Eigen::EigenSolver <Matrix_xd> eigenSolverMatrixValues(_matrixValues);

    auto stopTimeFindEigens = std::chrono::system_clock::now();
    auto durationTimeFindEigens = std::chrono::duration_cast<std::chrono::milliseconds>(
            stopTimeFindEigens - startTimeFindEigens);
    std::cout << "Duration to find eigenvectors and values: " << durationTimeFindEigens.count() << std::endl;

    _matrixValues.resize(0, 0); // Removes large matrix from memory; leads to faster write times for very large matrices

    auto startTimeSaveData = std::chrono::system_clock::now();

    MatrixOperations.save_data(_fileLocation, "eigenvectors_" + _fileNameNumSpins + ".csv",
                               eigenSolverMatrixValues.eigenvectors().real());
    MatrixOperations.save_data(_fileLocation, "eigenvalues_" + _fileNameNumSpins + ".csv",
                               eigenSolverMatrixValues.eigenvalues().imag());

    auto stopTimeSaveData = std::chrono::system_clock::now();
    auto durationTimeSaveData = std::chrono::duration_cast<std::chrono::milliseconds>(
            stopTimeSaveData - startTimeSaveData);
    std::cout << "Time to write to files: " << durationTimeSaveData.count() << std::endl;

    std::time_t stopTimeSaveData_cTimeUse = std::chrono::system_clock::to_time_t(stopTimeSaveData);
    std::cout << "\nFinished computation at: " << std::ctime(&stopTimeSaveData_cTimeUse) << std::endl;

    //std::cout << "Computing V * D * V^(-1) gives: " << std::endl << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

}