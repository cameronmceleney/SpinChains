//
// Created by Cameron Aidan McEleney on 20/11/2023.
//

#ifndef SPINCHAINS_IEFFECTIVEFIELD_H
#define SPINCHAINS_IEFFECTIVEFIELD_H

class IEffectiveField {
public:
    virtual double              EffectiveFieldX (const int& site, const int& layer, const double& mxLHS, const double& mxMID,
                                         const double& mxRHS, const double& dipoleTerm,
                                         const double& demagTerm, const double& current_time) = 0;

    // Description missing
    virtual double              EffectiveFieldY (const int& site, const int& layer, const double& myLHS, const double& myMID,
                                         const double& myRHS, const double& dipoleTerm,
                                         const double& demagTerm) = 0;

    // Description missing
    virtual double              EffectiveFieldZ (const int& site, const int& layer, const double& mzLHS, const double& mzMID,
                                         const double& mzRHS, const double& dipoleTerm,
                                         const double& demagTerm) = 0;

    ~IEffectiveField() = default;
};

#endif //SPINCHAINS_IEFFECTIVEFIELD_H
