#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"
#include "CommonLibs.h"

int main() {

    // GitHub Token: ***REMOVED*** (works as of 04 Jun 23)
    SpinChainEigenSolverClass SolverClass{};
    Numerical_Methods_Class NumericalMethods{};

    // Global file-related parameters
    GV.SetCurrentTime();
    GV.SetEmailWhenCompleted(false);

    // Global simulation parameters
    GV.SetAnisotropyField(0);
    GV.SetStaticBiasField(0.1);
    GV.SetNumSpins(3952);
    GV.SetExchangeMinVal(132.5);
    GV.SetExchangeMaxVal(132.5);
    GV.SetGyromagneticConstant(28.8);

    // Additional parameters and flags
    GV.SetIsFerromagnetic(true);
    GV.SetShouldFindEigenvalues(false);
    GV.SetIsExchangeUniform();
    std::string method = "RK2c";

    std::string outputFileID;
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> outputFileID;
    GV.SetFileNameBase("T"+outputFileID);

    GV.SetFilePath("windows");

    // I keep forgetting to check the exchanges, hence this warning
    if (GV.GetIsExchangeUniform())
        std::cout << "Uniform Exchange" << std::endl;
    else
        std::cout << "Non-Uniform Exchange" << std::endl;

    // Run the simulation
    if (GV.GetShouldFindEigenvalues()) {
        std::cout << "Finding eigenvalues and eigenvectors" << std::endl;
        SolverClass.CalculateEigenfrequencies(false);
    } else {
        NumericalMethods.NumericalMethodsMain();
        if (method == "RK2")
            NumericalMethods.SolveRK2();
        else if (method == "RK2c")
            NumericalMethods.SolveRK2Classic();
    }
    return 0;
}
