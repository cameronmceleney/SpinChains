#include <chrono>
#include <cmath>
#include <ctime>
#include "matrix_operations.h"
#include "common.h"

int main() {
    // TODO: rename variables
    // rename all variables following https://manual.gromacs.org/documentation/5.1-current/dev-manual/naming.html
    matrix_operations matrix_ops;

    int nspins; // number of spins in the chain
    double JI_min, JI_max; // minimum and maximum exchange (J) integral (I) values
    std::string file_location = R"(/Users/cameronmceleney/OneDrive - University of Glasgow/University/PhD/1st Year/C++/Chainspins/Data Outputs)";

    std::cout << "Enter the number of spins in the chain: ";
    std::cin >> nspins; // Takes user input for the number of spins

    /* Total number of spins (2N) is twice the number of spins (N) as there are two coupled equation (dx and dy)
     for each spin in the chain. */
    int total_equations = nspins * 2;

    // Creates unique filename by combining the number of spins with the keyword 'spins'
    std::string nspins_string;
    nspins_string = std::to_string(nspins) + "spins";

    std::cout << "Enter the Jmin value: ";
    std::cin >> JI_min;

    std::cout << "Enter the Jmax value: ";
    std::cin >> JI_max;

    // Compares user's exchange integrals and appends an explanatory string to the filename
    if (JI_min == JI_max){
        nspins_string += "-lin"; // Equal JI values means the system has linear exchange ('lin')
    }
    else if (JI_min != JI_max){
        nspins_string += "-nonlin"; // Unequal JI values means the system has nonlinear exchange ('nonlin')
    }
    else {
        /* no statement */
    }

    std::cout << "Filename is: " << nspins_string << std::endl; // Informs user of filename to enable directory searching via file explorer search function
    std::cout << "\nFiles will be outputted to: " << file_location << std::endl; // Showing the selected path will make it easier to find the file

    auto findeigs_start = std::chrono::system_clock::now(); // Separate start time (from system) for the computation of the eigenvalues
    std::time_t start_time = std::chrono::system_clock::to_time_t(findeigs_start);
    std::cout << "\nBegan computation at: "<< std::ctime(&start_time) << std::endl; // Useful for long computations where the start time may be forgotten

    MatrixcXd MatrixValues(total_equations,total_equations); // Generates the matrix but does not allocate memory. That is done as each element is calculated
    MatrixValues = matrix_ops.PopulateMatrix(nspins, JI_min, JI_max, 0.1, 29.2*2*M_PI);
    Eigen::EigenSolver<MatrixcXd> es(MatrixValues); //es stands for EigenSolver

    auto findeigs_stop = std::chrono::system_clock::now();
    auto findeigs_duration = std::chrono::duration_cast<std::chrono::milliseconds>(findeigs_stop - findeigs_start);
    std::cout << "Duration to find eigenvectors and values: " << findeigs_duration.count() << std::endl;

    MatrixValues.resize(0,0); // Removes large matrix from memory; leads to faster write times for very large matrices
    //std::cout << "Computing V * D * V^(-1) gives: " << std::endl << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

    auto savedata_start = std::chrono::system_clock::now();

    matrix_ops.SaveData(file_location,"eigenvectors_"+nspins_string+".csv", es.eigenvectors().real());
    matrix_ops.SaveData(file_location,"eigenvalues_"+nspins_string+".csv", es.eigenvalues().imag());

    auto savedata_stop = std::chrono::system_clock::now();
    auto savedata_duration = std::chrono::duration_cast<std::chrono::milliseconds>(savedata_stop - savedata_start);
    std::cout << "Time to write to files: " << savedata_duration.count() << std::endl;

    std::time_t end_time = std::chrono::system_clock::to_time_t(savedata_stop);
    std::cout << "\nFinished computation at: "<< std::ctime(&end_time) << std::endl;

    return 0;
}