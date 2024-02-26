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


    std::string baseFileName = "T" + simManager->outputFileID;
    int numString = 1;
    std::string letterString = "a";

    initialisationInstance->performInitialisation();
    configurationInstance->performInitialisation();

    simProgressBar.set_niter(168);
    simProgressBar.show_bar(true);

    for (double i = 11; i < 60; i+=0.2) {
        simProgressBar.update();
        // Add parameters to be changed here
        simParams->drivingFreq = i * 1e9;

        // For each new parameters set generate a new filename
        _generateNextFilename(baseFileName, letterString, numString);

        // Perform the simulation with new parameters
        initialisationInstance->reinitialise();
        configurationInstance->reinitialise();
        methodsInstance->performInitialisation();

        if ( i > 30) { i += 0.3; }
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

void SolversManager::_generateNextFilename(const std::string& baseName, std::string& letterSuffix, int& numSuffix) {
    if (simManager->hasNumericSuffix) {
        GV.SetFileNameBase(baseName + "_" + std::to_string(numSuffix++));
        return;
    } else {
        GV.SetFileNameBase(baseName + "_" + letterSuffix);

        // Logic for incrementing letterString.
        int j = letterSuffix.size() - 1;
        while ( j >= 0 ) {
            if ( letterSuffix[j] == 'z' ) {
                letterSuffix[j] = 'a';
                j--;
            } else {
                letterSuffix[j]++;
                break;
            }
        }

        if ( j < 0 )
            letterSuffix = 'a' + letterSuffix;
    }
}