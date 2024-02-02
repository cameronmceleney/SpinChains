//
// Created by Cameron McEleney on 02/02/2024.
//

// Corresponding header
#include "../include/SolversManager.h"

SolversManager::SolversManager(std::shared_ptr<SimulationManager> sharedSimManager,
                               std::shared_ptr<SimulationParameters> sharedSimParams,
                               std::shared_ptr<SimulationStates> sharedSimStates,
                               std::shared_ptr<SimulationFlags> sharedSimFlags)

    : SolversSuperClass(std::move(sharedSimManager), std::move(sharedSimParams), std::move(sharedSimStates), std::move(sharedSimFlags)) {
}

void SolversManager::Manager() {
    // TODO. Add method that automatically calls components which contain preset values (e.g. for maced2022breaking Fig.2)
    // Keep this set order to avoid errors

    if (simManager->massProduce) {
        massRunSimulations();
    } else
        singleSimulation();
}

void SolversManager::massRunSimulations() {

    auto initialisationInstance = SolversSuperClass::createSimulationInstance(simManager, simParams, simStates, simFlags);
    auto configurationInstance = SolversSuperClass::createConfigurationInstance(simManager, simParams, simStates, simFlags);
    auto methodsInstance = SolversSuperClass::createMethodsInstance(simManager, simParams, simStates, simFlags);


    progressbar mainBar(2, true);
    std::string baseFileName = "T" + simManager->outputFileID;
    int numString = 1;
    std::string letterString = "a";

    initialisationInstance->performInitialisation();
    configurationInstance->performInitialisation();

    simFlags->resetSimState = true;  // TODO must have this after initial configuration to ensure driving region etc is created


    for (int i = 2; i < 4; i++) {
        mainBar.update();
        simParams->drivingRegionWidth = i;

        if (simManager->hasNumericSuffix) {
            GV.SetFileNameBase(baseFileName + "_" + std::to_string(numString++));
        } else {
            GV.SetFileNameBase(baseFileName + "_" + letterString);

            // Logic for incrementing letterString.
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
                letterString = 'a' + letterString;
            }
        }

        configurationInstance->performInitialisation();
        methodsInstance->performInitialisation();

        // Assuming configurationInstance has a similar static function to create an instance and call performInitialisation
    }
}

void SolversManager::singleSimulation() {
    auto initialisationInstance = SolversSuperClass::createSimulationInstance(simManager, simParams, simStates, simFlags);
    auto configurationInstance = SolversSuperClass::createConfigurationInstance(simManager, simParams, simStates, simFlags);
    auto methodsInstance = SolversSuperClass::createMethodsInstance(simManager, simParams, simStates, simFlags);

    initialisationInstance->performInitialisation();
    configurationInstance->performInitialisation();
    methodsInstance->performInitialisation();
}