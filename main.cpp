#include "common.h"
#include "SpinChainEigenSolverClass.h"

int main() {

    int numberSpins; // number of spins in the chain

    std::cout << "Enter the number of spins in the chain: ";
    std::cin >> numberSpins; // Takes user input for the number of spins

    SpinChainEigenSolverClass SolverClass;
    void *output = SolverClass.DefineUserInputs(numberSpins);
    return 0;

    //std::cout << "Computing V * D * V^(-1) gives: " << std::endl << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

}