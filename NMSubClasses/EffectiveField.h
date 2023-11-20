//
// Created by Cameron McEleney on 31/10/2023.
//
#pragma once

#ifndef SPINCHAINS_EFFECTIVEFIELD_H
#define SPINCHAINS_EFFECTIVEFIELD_H

#include "NMMethods.h"
#include "IEffectiveField.h"

class EffectiveField : public NMMethods, public IEffectiveField {
public:
    EffectiveField(std::shared_ptr<SystemDataContainer> data) : NMMethods(data) {}
    ~EffectiveField() = default;
public:
    double              EffectiveFieldX (const int& site, const int& layer, const double& mxLHS, const double& mxMID,
                                         const double& mxRHS, const double& dipoleTerm,
                                         const double& demagTerm, const double& current_time) override;

    // Description missing
    double              EffectiveFieldY (const int& site, const int& layer, const double& myLHS, const double& myMID,
                                         const double& myRHS, const double& dipoleTerm,
                                         const double& demagTerm) override;

    // Description missing
    double              EffectiveFieldZ (const int& site, const int& layer, const double& mzLHS, const double& mzMID,
                                         const double& mzRHS, const double& dipoleTerm,
                                         const double& demagTerm) override;
};


#endif //SPINCHAINS_EFFECTIVEFIELD_H
