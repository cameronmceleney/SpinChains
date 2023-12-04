// C++ Standard Libraries
#include <iostream>

// C++ User Libraries (General)
#include "../include/CommonLibs.h"
#include "../include/SolversSuperClass.h"
// #include "../Other/working_on/SpinChainEigenSolverClass.h"

int main() {
    // GitHub Token: ***REMOVED*** (works as of 04 Jun 23)
    auto sharedSimParams = std::make_shared<SimulationParameters>();
    auto sharedSimStates = std::make_shared<SimulationStates>();
    auto sharedSimFlags = std::make_shared<SimulationFlags>();

    //SpinChainEigenSolverClass SolverClass{}; Needs to be turned into a separate class structure

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
    GV.SetDMIConstant(0);

    // Additional parameters and flags
    GV.SetIsFerromagnetic(true);
    GV.SetShouldFindEigenvalues(false);
    GV.SetIsExchangeUniform();

    std::string outputFileID;
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> outputFileID;
    GV.SetFileNameBase("T" + outputFileID);

    GV.SetFilePath("macos");

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
        auto initialisationInstance = SolversSuperClass::createSimulationInstance(sharedSimParams, sharedSimStates,
                                                                                  sharedSimFlags);
        auto configurationInstance = SolversSuperClass::createConfigurationInstance(sharedSimParams, sharedSimStates,
                                                                                    sharedSimFlags);
        auto methodsInstance = SolversSuperClass::createMethodsInstance(sharedSimParams, sharedSimStates,
                                                                        sharedSimFlags);

        initialisationInstance->performInitialisation();
        configurationInstance->performInitialisation();
        methodsInstance->performInitialisation();


    }
    return 0;
}
