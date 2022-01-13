#include "matrix_operations.h"

#include "linspace.h"

void matrix_operations::SaveData(std::string filePath, std::string fileName, Eigen::MatrixXd matrix)
{

    filePath_in = std::move(filePath);
    fileName_in = std::move(fileName);
    matrix_in = std::move(matrix);

    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(filePath_in+"/"+fileName_in);
    if (file.is_open())
    {
        file << matrix_in.format(CSVFormat);
        file.close();
    }
}

MatrixcXd matrix_operations::PopulateMatrix(int NumSpins, double ExchangeIntegral_min, double ExchangeIntegral_max, double BiasField, double GyroMagConst)
{

    NumSpins_in = NumSpins;
    ExchangeIntegral_min_in = ExchangeIntegral_min;
    ExchangeIntegral_max_in = ExchangeIntegral_max;
    BiasField_in = BiasField;
    GyroMagConst_in = GyroMagConst;

    auto const spins = static_cast<int>(NumSpins_in);
    auto const J_min = static_cast<double>(ExchangeIntegral_min_in);
    auto const J_max = static_cast<double>(ExchangeIntegral_max_in);
    auto const H0 = static_cast<double>(BiasField_in);
    auto const gamma = static_cast<double>(GyroMagConst_in);

    linspace Linspace{};
    std::vector<double> JI_values_linspace;
    Linspace.set_values(J_min, J_max, spins, true);
    JI_values_linspace = Linspace.generate_array(); //calls a linspace function similar to np.linspace()

    static double w = 0;
    int J_count = 0;
    //spinspairs = spins - 1
    const long max_N = spins * 2;

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
