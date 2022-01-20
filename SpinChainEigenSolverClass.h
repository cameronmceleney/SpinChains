#ifndef SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H
#define SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H

#include "CommonLibs.h"

class SpinChainEigenSolverClass {

private:
//  Dtype               Member Name                             //Comment
    double              _biasField;
    std::vector<double> _chainExchangeValues;
    std::string         _fileNameEigenSolver = GV.GetFileNameBase() + "spins"; // Creates unique filename by combining the number of spins with the keyword 'spins'
    Eigen::MatrixXd     _generatedMatrix;
    double              _gyroscopicMagneticConstant;
    Matrix_xd           _matrixValues;
    int                 _totalEquations;                        //Total number of spins (2*N) is twice the number of spins (N) as there are two coupled equation (dx and dy) for each spin in the chain.

public:
//  Dtype               Member Name                             //Comment
    void                CalculateEigFreqs();
    Matrix_xd           populate_matrix(double biasField, double gyroscopicMagneticConstant);
    void                save_data(std::string fileName, Matrix_xd generatedMatrix );
    void PrintVector(std::vector<double> inputVector);          // Writes a vector to the terminal
};

#endif //SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H
