#include "../include/SpinChainEigenSolverClass.h"
#include "../include/NMSuperClassTest.h"
#include "CommonLibs.h"

int main() {

    // GitHub Token: ***REMOVED*** (works as of 04 Jun 23)
    SpinChainEigenSolverClass SolverClass{};
    NMSuperClassTest NumericalMethods{};

    // Global file-related parameters
    GV.SetCurrentTime();
    GV.SetEmailWhenCompleted(false);

    // Global simulation parameters
    GV.SetAnisotropyField(0);
    GV.SetStaticBiasField(0.1);
    GV.SetNumSpins(4000);
    GV.SetExchangeMinVal(43.5);
    GV.SetExchangeMaxVal(132.0);
    GV.SetGyromagneticConstant(29.2);

    // Additional parameters and flags
    GV.SetIsFerromagnetic(true);
    GV.SetShouldFindEigenvalues(false);
    GV.SetIsExchangeUniform();

    std::string outputFileID;
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> outputFileID;
    GV.SetFileNameBase("T"+outputFileID);

    GV.SetFilePath("macos");

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
            std::shared_ptr<SharedVariableHolder> sharedVariables = std::make_shared<SharedVariableHolder>();
            NMSuperClassTest baseObj(sharedVariables);
            baseObj.executeDerivedMethod(); // Outputs: "Child method called!"
    }
    return 0;
}
