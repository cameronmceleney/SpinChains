//
// Created by Cameron McEleney on 31/10/2023.
//
#pragma once

#ifndef SPINCHAINS_EFFECTIVEFIELD_H
#define SPINCHAINS_EFFECTIVEFIELD_H

// Include the SimulationParameters
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"
#include "CommonLibs.h"

class EffectiveField {
private:
    SimulationParameters* _simParams; // Non-owning pointer to SimulationParameters
    SimulationStates* _simStates;
    SimulationFlags* _simFlags;
public:
    explicit EffectiveField(SimulationParameters* sharedSimParams,
                            SimulationStates* sharedSimStates,
                            SimulationFlags* sharedSimFlags);
    ~EffectiveField() = default;
    
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


#endif //SPINCHAINS_EFFECTIVEFIELD_H
