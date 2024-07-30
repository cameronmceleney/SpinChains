//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSINITIALISATION_H
#define SPINCHAINS_SOLVERSINITIALISATION_H

// C++ Standard Libraries
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

// C++ User Libraries (General)
#include "GlobalVariables.h"

// C++ User Libraries (Class' Parent)
#include "SolversSuperClass.h"

class SolversInitialisation: public SolversSuperClass {
private:
    const double        _BOHR_MAGNETON = 9.274e-24;

    void                _setSimulationFlags();
    void                _setSimulationParameters();
    void                _generateRemainingParameters();
    void                _setMaterialParameters();
    void                _guardClauses();

    void                Initialise();
    void                _recalculateParameters();

public:
    SolversInitialisation(std::shared_ptr<SimulationManager> sharedSimManager,
                          std::shared_ptr<SimulationParameters> sharedSimParams,
                     std::shared_ptr<SimulationStates> sharedSimStates,
                     std::shared_ptr<SimulationFlags> sharedSimFlags);
    ~SolversInitialisation() = default;
public:
    void                testModifyingDouble(double  newValue);
    void                performInitialisation() override { Initialise(); };
    void                reinitialise() override { _recalculateParameters(); };
};


#endif //SPINCHAINS_SOLVERSINITIALISATION_H
