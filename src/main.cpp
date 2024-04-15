
// C++ Standard Libraries


// C++ User Libraries (General)
//#include "../checks_examples/check_for_metal_device.h"
//#include "../test/mappingSpeedTests.hpp"
#include "../libs/CommonDefinitions.h"
#include "../include/SolversSuperClass.h"
// #include "../Other/working_on/SpinChainEigenSolverClass.h"

int main() {
    //CoffeeExample coffeeExample{};
    //const int numElements = 1 << 18;
    //const bool debugMode = false;
    //std::string methodType = "simpleVectorSineAdd";
    //coffeeExample.vectorAddition(numElements, methodType, debugMode);

    auto sharedSimManager = std::make_shared<SimulationManager>();

    //SpinChainEigenSolverClass SolverClass{}; Needs to be turned into a separate class structure

    // Global file-related parameters
    GV.SetCurrentTime();
    GV.SetEmailWhenCompleted(false);

    // Global simulation parameters
    GV.SetAnisotropyField(0);
    GV.SetStaticBiasField(0.1);
    GV.SetNumSpins(1400);
    GV.SetExchangeMinVal(43.5);
    GV.SetExchangeMaxVal(132);
    GV.SetGyromagneticConstant(29.2);
    GV.SetDMIConstant(0);  // use negative to flip to match python for now

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

        sharedSimManager->massProduce = false;
        sharedSimManager->hasNumericSuffix = true;

        auto managerInstance = SolversSuperClass::createSimulationManager(sharedSimManager, sharedSimParams, sharedSimStates,
                                                                          sharedSimFlags);
        managerInstance->performInitialisation();
    }
    return 0;
}
