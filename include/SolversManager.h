//
// Created by Cameron McEleney on 02/02/2024.
//

#ifndef SPINCHAINS_SOLVERSMANAGER_H
#define SPINCHAINS_SOLVERSMANAGER_H

// C++ Standard Libraries


// C++ User Libraries (General)
#include "GlobalVariables.h"

// C++ User Libraries (Class' Parent)
#include "SolversSuperClass.h"

class SolversManager: public SolversSuperClass {
private:
    void                Manager();
    void                massRunSimulations();
    void                singleSimulation();
    void                _generateNextFilename(const std::string& baseName, std::string& letterSuffix, int& numSuffix);

public:
    SolversManager( std::shared_ptr<SimulationManager> sharedSimManager,
                    std::shared_ptr<SimulationParameters> sharedSimParams,
                    std::shared_ptr<SimulationStates> sharedSimStates,
                    std::shared_ptr<SimulationFlags> sharedSimFlags);
    ~SolversManager() = default;
public:
    void                performInitialisation() override { Manager(); };
    void                reinitialise() override { Manager(); };
};


#endif //SPINCHAINS_SOLVERSMANAGER_H
