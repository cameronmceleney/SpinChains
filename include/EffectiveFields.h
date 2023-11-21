//
// Created by Cameron McEleney on 31/10/2023.
//
#pragma once

#ifndef SPINCHAINS_EFFECTIVEFIELDS_H
#define SPINCHAINS_EFFECTIVEFIELDS_H

// C++ Standard Library
#include <vector>

// C++ User Libraries (General)
#include "CommonLibs.h"

// C++ User Libraries (General)
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"

class EffectiveFields {
private:
    SimulationParameters* _simParams; // Non-owning pointer to SimulationParameters
    SimulationStates* _simStates;
    SimulationFlags* _simFlags;
public:
    explicit EffectiveFields(SimulationParameters* sharedSimParams,
                            SimulationStates* sharedSimStates,
                            SimulationFlags* sharedSimFlags);
    ~EffectiveFields() = default;
    
public:
    double              EffectiveFieldX (const int& site, const int& layer, const double& mxLHS, const double& mxMID,
                                         const double& mxRHS, const double& dipoleTerm,
                                         const double& demagTerm, const double& current_time);

    // Description missing
    double              EffectiveFieldY (const int& site, const int& layer, const double& myLHS, const double& myMID,
                                         const double& myRHS, const double& dipoleTerm,
                                         const double& demagTerm);

    // Description missing
    double              EffectiveFieldZ (const int& site, const int& layer, const double& mzLHS, const double& mzMID,
                                         const double& mzRHS, const double& dipoleTerm,
                                         const double& demagTerm);
};


#endif //SPINCHAINS_EFFECTIVEFIELDS_H
