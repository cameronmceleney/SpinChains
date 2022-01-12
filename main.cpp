#include <chrono>
#include <cmath>
#include <complex>
#include <ctime>
#include <iostream>
#include <sstream>

#include "linspace.h"
#include "saveMatrix.h"
#include "common.h"





template<typename S>
MatrixcXd PopulateMatrix(int num_spins_in_chain, S ExchangeIntegral_min, S ExchangeIntegral_max, S BiasField, S GyroMagConst) {

    //std::vector<double> JI_values_linspace = linspace(J_min, JI_max, num_spinpairs); //function is like np.linspace()

    auto const spins = static_cast<long>(num_spins_in_chain);
    auto const J_min = static_cast<double>(ExchangeIntegral_min);
    auto const J_max = static_cast<double>(ExchangeIntegral_max);
    auto const H0 = static_cast<double>(BiasField);
    auto const gamma = static_cast<double>(GyroMagConst);

    linspace Linspace;
    std::vector<double> JI_values_linspace;
    Linspace.set_values(J_min, J_max, spins, true);
    JI_values_linspace = Linspace.generate_array();

    static double w = 0;
    int J_count = 0;
    const long spinpairs = spins - 1, max_N = spins * 2;
    //std::vector<double> JI_values_linspace = myObj.findlinspace(J_min, J_max, spinpairs); //calls a linspace function similar to np.linspace()
    //adds a zero to the start and end of JI_values. Insert is faster for large values of numbers compared to push_back()
    std::vector<double> J_array{0}; //initialised with a zero to account for the (P-1)th spin
    J_array.insert(J_array.end(), JI_values_linspace.begin(), JI_values_linspace.end());
    J_array.push_back(0); //add a zero to the end to account for the (N+1)th spin

    MatrixcXd matrix(max_N, max_N);
    matrix.setZero();

    for (long row_num = 0; row_num < max_N; row_num++) {

        if (row_num % 2 == 0) {
            //These are even numbered rows, containing the dx/dt terms, so the 'Odd Rule' applies (see notes)

            if (row_num == 0) {
                //exception for the first dx/dt row (first row overall), as there is no spin on the LHS of this position
                //static_cast<long double>(-1.0)
                matrix(row_num,0) = -1.0  * w / gamma;
                matrix(row_num,1) = (J_array[J_count] + J_array[J_count + 1] + H0);
                matrix(row_num,3) = -1.0 * J_array[J_count + 1];

            } else if (row_num > 0 and row_num < max_N - 2) {
                //handles all other even numbered rows
                matrix(row_num,row_num - 1) = -1.0 * J_array[J_count];
                matrix(row_num,row_num + 0) = -1.0  * w / gamma;
                matrix(row_num,row_num + 1) = (J_array[J_count] + J_array[J_count + 1] + H0);
                matrix(row_num,row_num + 3) = -1.0 * J_array[J_count + 1];

            } else if (row_num == max_N - 2) {
                // exception for the final dx/dt row (penultimate row overall), as there is no spin on the RHS of this position

                matrix(row_num,max_N - 1) = (J_array[J_count] + J_array[J_count + 1] + H0);
                matrix(row_num,max_N - 2) = -1.0  * w / gamma;
                matrix(row_num,max_N - 3) = -1.0 * J_array[J_count];

            } else {

                std::cout << "Error with generating the dx/dt terms on row #{row_num}. Exiting..." << std::endl;
                std::exit(1);

            }
        }

        else if (row_num % 2 == 1) {
            //These are odd numbered rows, containing the dy/dt terms, so the 'Even Rule' applies (see notes)
            if (row_num == 1) {
                //exception for the first dy/dt row (second row overall), as there is no spin on the LHS of this position
                matrix(row_num,0) = -1.0 * (J_array[J_count] + J_array[J_count + 1] + H0);
                matrix(row_num,1) = -1.0  * w / gamma;
                matrix(row_num,2) = J_array[J_count + 1];

            } else if (row_num > 1 and row_num < max_N - 1) {
                //handles all other odd numbered rows
                matrix(row_num,row_num - 3) = J_array[J_count];
                matrix(row_num,row_num - 1) = -1.0 * (J_array[J_count] + J_array[J_count + 1] + H0);
                matrix(row_num,row_num + 0) = -1.0  * w / gamma;
                matrix(row_num,row_num + 1) = J_array[J_count + 1];
            } else if (row_num == max_N - 1) {
                //exception for the final dy/dt row (final row overall), as there is no spin on the RHS of this position
                matrix(row_num,max_N - 1) = -1.0  * w / gamma;
                matrix(row_num,max_N - 2) = -1.0 * (J_array[J_count] + J_array[J_count + 1] + H0);
                matrix(row_num,max_N - 4) = J_array[J_count];
            } else {

                std::cout << "Error with generating the dy/dt terms on row #{row_num}. Exiting..." << std::endl;
                std::exit(1);
            }
            J_count += 1;
        }

        else {
            std::cout << "An issue arose in generating the matrix. Exiting..." << std::endl;
            std::exit(1);
        }
    }

    //std::cout << matrix << "\n\n" << std::endl;
    matrix *= gamma;
    return matrix;
}

int main() {

    saveMatrix savematrix;
    long nspins;
    double JI_min, JI_max;
    std::string file_location = "C:\\Users\\Cameron McEleney\\OneDrive - University of Glasgow\\University\\PhD\\1st Year\\C++\\Chainspins\\Data Outputs";

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

    auto findeigs_start = std::chrono::high_resolution_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(findeigs_start);
    std::cout << "\nBegan computation at: "<< std::ctime(&start_time) << std::endl;

    MatrixcXd MatrixValues(total_N,total_N);
    MatrixValues = PopulateMatrix(nspins, JI_min, JI_max, 0.1, 29.2*2*M_PI);
    Eigen::EigenSolver<MatrixcXd> es(MatrixValues); //ces stands for ComplexEigenSolver for ease

    auto findeigs_stop = std::chrono::high_resolution_clock::now();
    auto findeigs_duration = std::chrono::duration_cast<std::chrono::milliseconds>(findeigs_stop - findeigs_start);
    std::cout << "Duration to find eigenvectors and values: " << findeigs_duration.count() << std::endl;

    MatrixValues.resize(0,0); //to remove from memory by resizing to 0
    //std::cout << "Computing V * D * V^(-1) gives: " << std::endl << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

    auto savedata_start = std::chrono::high_resolution_clock::now();

    savematrix.saveData(file_location,"eigenvectors_"+nspins_string+".csv", es.eigenvectors().real());
    savematrix.saveData(file_location,"eigenvalues_"+nspins_string+".csv", es.eigenvalues().imag());

    auto savedata_stop = std::chrono::high_resolution_clock::now();
    auto savedata_duration = std::chrono::duration_cast<std::chrono::milliseconds>(savedata_stop - savedata_start);
    std::cout << "Time to write to files: " << savedata_duration.count() << std::endl;

    std::time_t end_time = std::chrono::system_clock::to_time_t(savedata_stop);
    std::cout << "\nFinished computation at: "<< std::ctime(&end_time) << std::endl;

    return 0;
}