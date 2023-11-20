//
// Created by Cameron McEleney on 31/10/2023.
//

#include "../include/EffectiveField.h"

EffectiveField::EffectiveField(SystemDataContainer* data) : systemData(data) {}


double EffectiveField::EffectiveFieldX(const int& site, const int& layer, const double& mxLHS, const double& mxMID,
                                                const double& mxRHS, const double& dipoleTerm,
                                                const double& demagTerm, const double& current_time) {
    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hx;

    if (systemData->isFm) {
        if (site >= systemData->drivingRegionLhs && site <= systemData->drivingRegionRhs) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (systemData->driveAllLayers || layer == 0)
                hx = systemData->exchangeVec[site - 1] * mxLHS + systemData->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm
                        + systemData->dynamicBiasField * cos(systemData->drivingAngFreq * current_time);
            else if  (systemData->hasStaticDrive)
                hx = systemData->exchangeVec[site - 1] * mxLHS + systemData->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm
                        +systemData->dynamicBiasField;
            else if  ((!systemData->driveAllLayers && layer != 0))
                hx = systemData->exchangeVec[site - 1] * mxLHS + systemData->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm;
        } else
            // All spins along x which are not within the driving region
            hx = systemData->exchangeVec[site - 1] * mxLHS + systemData->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm;
    } else if (!systemData->isFm) {
        if (site >= systemData->drivingRegionLhs && site <= systemData->drivingRegionRhs) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (systemData->hasStaticDrive)
                hx = -1.0 * (systemData->exchangeVec[site - 1] * mxLHS + systemData->exchangeVec[site] * mxRHS + systemData->dynamicBiasField);
            else if (!systemData->hasStaticDrive)
                hx = -1.0 * (systemData->exchangeVec[site - 1] * mxLHS + systemData->exchangeVec[site] * mxRHS) + systemData->dynamicBiasField * cos(systemData->drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx = -1.0 * (systemData->exchangeVec[site - 1] * mxLHS + systemData->exchangeVec[site] * mxRHS);
    }

    return hx;
}
double EffectiveField::EffectiveFieldY(const int& site, const int& layer, const double& myLHS, const double& myMID, const double& myRHS,
                                                const double &dipoleTerm,  const double& demagTerm) {
    // The effective field (H_eff) y-component acting upon a given magnetic moment (site), abbreviated to 'hy'
    double hy;

    if (systemData->isFm) {
        hy = systemData->exchangeVec[site-1] * myLHS + systemData->exchangeVec[site] * myRHS + dipoleTerm + demagTerm;
    } else if (!systemData->isFm) {
        hy = -1.0 * (systemData->exchangeVec[site-1] * myLHS + systemData->exchangeVec[site] * myRHS);
    }

    return hy;
}
double EffectiveField::EffectiveFieldZ(const int& site, const int& layer, const double& mzLHS, const double& mzMID, const double& mzRHS,
                                                const double& dipoleTerm, const double& demagTerm) {
    // The effective field (H_eff) z-component acting upon a given magnetic moment (site), abbreviated to 'hz'
    double hz;

    if (systemData->isFm) {
        hz = systemData->exchangeVec[site-1] * mzLHS + systemData->exchangeVec[site] * mzRHS + dipoleTerm + demagTerm
                + GV.GetStaticBiasField();
    } else if (!systemData->isFm) {
        if (mzMID > 0)
            hz = GV.GetStaticBiasField() + systemData->anisotropyField - (systemData->exchangeVec[site-1] * mzLHS + systemData->exchangeVec[site] * mzRHS);
        else if (mzMID < 0)
            hz = GV.GetStaticBiasField() - systemData->anisotropyField - (systemData->exchangeVec[site-1] * mzLHS + systemData->exchangeVec[site] * mzRHS);
    }

    return hz;
}