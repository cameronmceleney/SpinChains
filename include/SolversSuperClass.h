//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSSUPERCLASS_H
#define SPINCHAINS_SOLVERSSUPERCLASS_H

// C++ Standard Libraries
//#include <algorithm>  // Will need in future. Left here so I don't forget

// C++ Third Party Libraries

// C++ User Libraries (General)
#include "../libs/CommonDefinitions.h"
#include "../libs/CommonStructures.h"
#include "../libs/progressbar.hpp"

// C++ User Libraries (Containers)
#include "SimulationManager.h"
#include "SimulationFlags.h"
#include "SimulationParameters.h"
#include "SimulationStates.h"

class SolversSuperClass {
protected:
    // Getter for child classes to access sharedVariables
    std::shared_ptr<SimulationManager> simManager;
    std::shared_ptr<SimulationParameters> simParams;
    std::shared_ptr<SimulationStates> simStates;
    std::shared_ptr<SimulationFlags> simFlags;
    CommonStructures::Timer methodTimer;
    progressbar simProgressBar;

public:
    //SolversSuperClass(std::shared_ptr<SimulationParameters> simParams);
    SolversSuperClass( std::shared_ptr<SimulationManager> sharedSimManager,
                       std::shared_ptr<SimulationParameters> sharedSimParams,
                       std::shared_ptr<SimulationStates> sharedSimStates,
                       std::shared_ptr<SimulationFlags> sharedSimFlags );

    virtual ~SolversSuperClass(); // Destructor (though not strictly necessary with smart pointers)

public:
    virtual void performInitialisation() = 0;
    virtual void reinitialise() = 0;

    static std::shared_ptr<SolversSuperClass> 
    createSimulationManager( std::shared_ptr<SimulationManager> sharedSimManager,
                             std::shared_ptr<SimulationParameters> sharedSimParams,
                             std::shared_ptr<SimulationStates> sharedSimStates,
                             std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SolversSuperClass>
    createSimulationInstance( std::shared_ptr<SimulationManager> sharedSimManager,
                              std::shared_ptr<SimulationParameters> sharedSimParams,
                              std::shared_ptr<SimulationStates> sharedSimStates,
                              std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SolversSuperClass>
    createConfigurationInstance( std::shared_ptr<SimulationManager> sharedSimManager,
                                 std::shared_ptr<SimulationParameters> sharedSimParams,
                                 std::shared_ptr<SimulationStates> sharedSimStates,
                                 std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SolversSuperClass>
    createDataHandlingInstance( std::shared_ptr<SimulationManager> sharedSimManager,
                                std::shared_ptr<SimulationParameters> sharedSimParams,
                                std::shared_ptr<SimulationStates> sharedSimStates,
                                std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SolversSuperClass>
    createMethodsInstance( std::shared_ptr<SimulationManager> sharedSimManager,
                           std::shared_ptr<SimulationParameters> sharedSimParams,
                           std::shared_ptr<SimulationStates> sharedSimStates,
                           std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SimulationManager> getSimulationManagerContainer();
    
    static std::shared_ptr<SimulationParameters> getSimulationParamsContainer();

    static std::shared_ptr<SimulationStates> getSimulationStatesContainer();

    static std::shared_ptr<SimulationFlags> getSimulationFlagsContainer();

    static std::pair<double, double> testInstance();
};


#endif //SPINCHAINS_SOLVERSSUPERCLASS_H
