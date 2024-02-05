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
        simFlags->resetSimState = true;
        massRunSimulations();
    } else
        singleSimulation();
}

void SolversManager::massRunSimulations() {

    auto initialisationInstance = SolversSuperClass::createSimulationInstance(simManager, simParams, simStates, simFlags);
    auto configurationInstance = SolversSuperClass::createConfigurationInstance(simManager, simParams, simStates, simFlags);
    auto methodsInstance = SolversSuperClass::createMethodsInstance(simManager, simParams, simStates, simFlags);


    progressbar mainBar(68, true);
    std::string baseFileName = "T" + simManager->outputFileID;
    int numString = 1;
    std::string letterString = "a";

    initialisationInstance->performInitialisation();
    configurationInstance->performInitialisation();

    for (int i = 132; i < 201; i++) {
        mainBar.update();
        // Add parameters to be changed here
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