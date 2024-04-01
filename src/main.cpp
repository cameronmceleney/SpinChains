// C++ Standard Libraries
#include <iostream>

// C++ User Libraries (General)
//#include "../test/mappingSpeedTests.hpp"
#include "../libs/CommonDefinitions.h"
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
    GV.SetStaticBiasField(0.005);
    GV.SetNumSpins(8200);
    GV.SetExchangeMinVal(52.0);
    GV.SetExchangeMaxVal(52.0);
    GV.SetGyromagneticConstant(29.2);
    GV.SetDMIConstant(1.0);  // use negative to flip to match python for now

    // Additional parameters and flags
    GV.SetIsFerromagnetic(true);
    GV.SetShouldFindEigenvalues(false);
    GV.SetIsExchangeUniform();

    std::string outputFileID;
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> outputFileID;
    GV.SetFileNameBase("T" + outputFileID);
    sharedSimManager->outputFileID = outputFileID;

    // OS Detection and setting the file path accordingly
    #if defined(_WIN32) || defined(_WIN64)
        GV.SetFilePath("windows");
    #elif defined(__linux__)
        GV.SetFilePath("linux");
    #elif defined(__APPLE__) && defined(__MACH__)
        GV.SetFilePath("macos");
    #else
        GV.SetFilePath("other");
    #endif

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

        sharedSimManager->massProduce = true;
        sharedSimManager->hasNumericSuffix = true;

        auto managerInstance = SolversSuperClass::createSimulationManager(sharedSimManager, sharedSimParams, sharedSimStates,
                                                                          sharedSimFlags);
        managerInstance->performInitialisation();
    }
    return 0;
}
