#include <chrono>
#include <cmath>
#include <ctime>

#include "matrix_operations.h"
#include "common.h"

int main() {

    matrix_operations matrix_ops;
    int nspins;
    double JI_min, JI_max;
    std::string file_location = R"(/Users/cameronmceleney/OneDrive - University of Glasgow/University/PhD/1st Year/C++/Chainspins/Data Outputs)";

    std::cout << "Enter the number of spins in the chain: ";
    std::cin >> nspins;

    int total_N = nspins * 2;
    std::string nspins_string;
    nspins_string = std::to_string(nspins) + "spins";

    std::cout << "Enter the Jmin value: ";
    std::cin >> JI_min;

    std::cout << "Enter the Jmax value: ";
    std::cin >> JI_max;

    std::cout << "File extension is: " << nspins_string << std::endl;
    std::cout << "\nFiles will be outputted to: " << file_location << std::endl;

    auto findeigs_start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(findeigs_start);
    std::cout << "\nBegan computation at: "<< std::ctime(&start_time) << std::endl;

    MatrixcXd MatrixValues(total_N,total_N);
    MatrixValues = matrix_ops.PopulateMatrix(nspins, JI_min, JI_max, 0.1, 29.2*2*M_PI);
    Eigen::EigenSolver<MatrixcXd> es(MatrixValues); //ces stands for ComplexEigenSolver for ease

    auto findeigs_stop = std::chrono::system_clock::now();
    auto findeigs_duration = std::chrono::duration_cast<std::chrono::milliseconds>(findeigs_stop - findeigs_start);
    std::cout << "Duration to find eigenvectors and values: " << findeigs_duration.count() << std::endl;

    MatrixValues.resize(0,0); //to remove from memory by resizing to 0
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