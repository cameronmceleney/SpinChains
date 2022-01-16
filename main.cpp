#include "common.h"
#include "SpinChainEigenSolverClass.h"

int main() {

    int numberSpins; // number of spins in the chain

    std::cout << "Enter the number of spins in the chain: ";
    std::cin >> numberSpins; // Takes user input for the number of spins

    SpinChainEigenSolverClass SolverClass;
    SolverClass.SolveInputs(numberSpins);

    return 0;

}