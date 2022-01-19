#ifndef SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H
#define SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H

#include "CommonLibs.h"

class SpinChainEigenSolverClass {

private:
    int _totalEquations; //Total number of spins (2*N) is twice the number of spins (N) as there are two coupled equation (dx and dy) for each spin in the chain.

    Eigen::MatrixXd _generatedMatrix;
    std::vector<double> _chainExchangeValues;

    double _biasField;
    double _gyroscopicMagneticConstant;

    std::string _fileNameEigenSolver = GV.GetFileNameBase() + "spins"; // Creates unique filename by combining the number of spins with the keyword 'spins'

    Matrix_xd _matrixValues;

public:

    void SolveInputs();
    Matrix_xd populate_matrix(double biasField, double gyroscopicMagneticConstant);
    void save_data(std::string fileName, Matrix_xd generatedMatrix );
    // Writes a vector to the terminal
    void PrintVector(std::vector<double> inputVector);
};

#endif //SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H
