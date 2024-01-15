// C++ Standard Libraries
#include <iostream>

// C++ User Libraries (General)
#include "../include/CommonLibs.h"
#include "../include/SolversSuperClass.h"
// #include "../Other/working_on/SpinChainEigenSolverClass.h"
#include "../libs/progressbar.hpp"

int main() {
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
    GV.SetExchangeMinVal(150);
    GV.SetExchangeMaxVal(150);
    GV.SetGyromagneticConstant(28);
    GV.SetDMIConstant(15);

    // Additional parameters and flags
    GV.SetIsFerromagnetic(true);
    GV.SetShouldFindEigenvalues(false);
    GV.SetIsExchangeUniform();

    std::string outputFileID;
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> outputFileID;
    GV.SetFileNameBase("T" + outputFileID);

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
        auto initialisationInstance = SolversSuperClass::createSimulationInstance(sharedSimParams, sharedSimStates,
                                                                                  sharedSimFlags);
        auto configurationInstance = SolversSuperClass::createConfigurationInstance(sharedSimParams, sharedSimStates,
                                                                                    sharedSimFlags);
        auto methodsInstance = SolversSuperClass::createMethodsInstance(sharedSimParams, sharedSimStates,
                                                                        sharedSimFlags);

        initialisationInstance->performInitialisation();
        configurationInstance->performInitialisation();
        sharedSimFlags->resetSimState = true;

        std::string letterString = "a";
        progressbar mainBar(50);
        for (int i = 1; i < 50; i++) {
            mainBar.update();
            sharedSimParams->drivingRegionWidth = i;

            GV.SetFileNameBase(outputFileID + letterString);

            // Increment LETTER
            int j = letterString.size() - 1;
            while (j >= 0) {
                if (letterString[j] == 'z') {
                    letterString[j] = 'a';
                    j--;
                } else {
                    letterString[j]++;
                    break;
                }
            }
            if (j < 0) {
                letterString = 'a' + letterString; // Prepend 'a' when all characters were 'z'
            }

            methodsInstance->performInitialisation();
            configurationInstance->performInitialisation();
        }
    }
    return 0;
}
