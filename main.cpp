#include "common.h"
#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"

int main() {

    int numberSpins; // number of spins in the chain

    std::cout << "Enter the number of spins in the chain: ";
    std::cin >> numberSpins; // Takes user input for the number of spins

    SpinChainEigenSolverClass SolverClass;
    Numerical_Methods_Class NMMethods;

    //SolverClass.SolveInputs(numberSpins);
    NMMethods.RK2(numberSpins);

    return 0;

}