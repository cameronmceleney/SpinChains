#include <chrono>
#include <cmath>
#include <ctime>
#include "matrix_operations.h"
#include "common.h"

int main() {
    // TODO: rename variables
    // rename all variables following https://manual.gromacs.org/documentation/5.1-current/dev-manual/naming.html
    matrix_operations MatrixOps;

    int g_numberSpins; // number of spins in the chain
    double exchangeMinimum, exchangeMaximum; // minimum and maximum exchange (J) integral (I) values
    std::string file_location = R"(/Users/cameronmceleney/OneDrive - University of Glasgow/University/PhD/1st Year/C++/Chainspins/Data Outputs)";

    std::cout << "Enter the number of spins in the chain: ";
    std::cin >> g_numberSpins; // Takes user input for the number of spins

    /* Total number of spins (2N) is twice the number of spins (N) as there are two coupled equation (dx and dy)
     for each spin in the chain. */
    int const c_totalEquations = g_numberSpins * 2;

    // Creates unique filename by combining the number of spins with the keyword 'spins'
    std::string fileNameNumSpins;
    fileNameNumSpins = std::to_string(g_numberSpins) + "spins";

    std::cout << "Enter the Jmin value: ";
    std::cin >> exchangeMinimum;

    std::cout << "Enter the Jmax value: ";
    std::cin >> exchangeMaximum;

    // Compares user's exchange integrals and appends an explanatory string to the filename
    if (exchangeMinimum == exchangeMaximum){
        fileNameNumSpins += "-lin"; // Equal JI values means the system has linear exchange ('lin')
    }
    else if (exchangeMinimum != exchangeMaximum){
        fileNameNumSpins += "-nonlin"; // Unequal JI values means the system has nonlinear exchange ('nonlin')
    }
    else {
        /* no statement */
    }

    std::cout << "Filename is: " << fileNameNumSpins << std::endl; // Informs user of filename to enable directory searching via file explorer search function
    std::cout << "\nFiles will be outputted to: " << file_location << std::endl; // Showing the selected path will make it easier to find the file

    auto startTimeFindEigens = std::chrono::system_clock::now(); // Separate start time (from system) for the computation of the eigenvalues
    std::time_t startTimeFindEigens_cTimeUse = std::chrono::system_clock::to_time_t(startTimeFindEigens);
    std::cout << "\nBegan computation at: "<< std::ctime(&startTimeFindEigens_cTimeUse) << std::endl; // Useful for long computations where the start time may be forgotten

    MatrixcXd MatrixValues(c_totalEquations,c_totalEquations); // Generates the matrix but does not allocate memory. That is done as each element is calculated
    MatrixValues = MatrixOps.PopulateMatrix(g_numberSpins, exchangeMinimum, exchangeMaximum, 0.1, 29.2*2*M_PI);
    Eigen::EigenSolver<MatrixcXd> eigenSolverMatrixValues(MatrixValues);

    auto stopTimeFindEigens = std::chrono::system_clock::now();
    auto durationTimeFindEigens = std::chrono::duration_cast<std::chrono::milliseconds>(stopTimeFindEigens - startTimeFindEigens);
    std::cout << "Duration to find eigenvectors and values: " << durationTimeFindEigens.count() << std::endl;

    MatrixValues.resize(0,0); // Removes large matrix from memory; leads to faster write times for very large matrices
    //std::cout << "Computing V * D * V^(-1) gives: " << std::endl << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

    auto startTimeSaveData = std::chrono::system_clock::now();

    MatrixOps.SaveData(file_location,"eigenvectors_"+fileNameNumSpins+".csv", eigenSolverMatrixValues.eigenvectors().real());
    MatrixOps.SaveData(file_location,"eigenvalues_"+fileNameNumSpins+".csv", eigenSolverMatrixValues.eigenvalues().imag());

    auto stopTimeSaveData = std::chrono::system_clock::now();
    auto durationTimeSaveData = std::chrono::duration_cast<std::chrono::milliseconds>(stopTimeSaveData - startTimeSaveData);
    std::cout << "Time to write to files: " << durationTimeSaveData.count() << std::endl;

    std::time_t stopTimeSaveData_cTimeUse = std::chrono::system_clock::to_time_t(stopTimeSaveData);
    std::cout << "\nFinished computation at: "<< std::ctime(&stopTimeSaveData_cTimeUse) << std::endl;

    return 0;
}