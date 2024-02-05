// C++ Standard Libraries
#include <iostream>

// C++ User Libraries (General)
#include "../include/CommonLibs.h"
#include "../include/SolversSuperClass.h"
// #include "../Other/working_on/SpinChainEigenSolverClass.h"

int main() {
    auto sharedSimManager = std::make_shared<SimulationManager>();

    //SpinChainEigenSolverClass SolverClass{}; Needs to be turned into a separate class structure

    // Global file-related parameters
    GV.SetCurrentTime();
    GV.SetEmailWhenCompleted(false);

    // Global simulation parameters
    GV.SetAnisotropyField(0);
    GV.SetStaticBiasField(0.1);
    GV.SetNumSpins(3400);
    GV.SetExchangeMinVal(32.5);
    GV.SetExchangeMaxVal(32.5);
    GV.SetGyromagneticConstant(28.0);
    GV.SetDMIConstant(2.5);

    // Additional parameters and flags
    GV.SetIsFerromagnetic(true);
    GV.SetShouldFindEigenvalues(false);
    GV.SetIsExchangeUniform();

    std::string outputFileID;
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> outputFileID;
    GV.SetFileNameBase("T" + outputFileID);
    sharedSimManager->outputFileID = outputFileID;

    GV.SetFilePath("windows");

    GV.SetNumericalMethod("RK2p");

    // I keep forgetting to check the exchanges, hence this warning
    if ( GV.GetIsExchangeUniform())
        std::cout << "Uniform Exchange" << std::endl;
    else
        std::cout << "Non-Uniform Exchange" << std::endl;

    // Run the simulation
    if ( GV.GetShouldFindEigenvalues()) {
        std::cout << "Finding eigenvalues and eigenvectors" << std::endl;
    } else {
        auto sharedSimParams = std::make_shared<SimulationParameters>();
        auto sharedSimStates = std::make_shared<SimulationStates>();
        auto sharedSimFlags = std::make_shared<SimulationFlags>();

        sharedSimManager->massProduce = false;
        sharedSimManager->hasNumericSuffix = true;

        auto managerInstance = SolversSuperClass::createSimulationManager(sharedSimManager, sharedSimParams, sharedSimStates,
                                                                          sharedSimFlags);
        managerInstance->performInitialisation();
    }
    return 0;
}
