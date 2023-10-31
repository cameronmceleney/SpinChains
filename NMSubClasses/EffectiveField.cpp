//
// Created by Cameron McEleney on 31/10/2023.
//

#include "EffectiveField.h"

double EffectiveField::EffectiveFieldX(const int& site, const int& layer, const double& mxLHS, const double& mxMID,
                                                const double& mxRHS, const double& dipoleTerm,
                                                const double& demagTerm, const double& current_time) {
    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hx;

    if (isFm) {
        if (site >= drivingRegionLhs && site <= drivingRegionRhs) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (driveAllLayers || layer == 0)
                hx = exchangeVec[site - 1] * mxLHS + exchangeVec[site] * mxRHS + dipoleTerm + demagTerm
                        + dynamicBiasField * cos(drivingAngFreq * current_time);
            else if  (hasStaticDrive)
                hx = exchangeVec[site - 1] * mxLHS + exchangeVec[site] * mxRHS + dipoleTerm + demagTerm
                        +dynamicBiasField;
            else if  ((!driveAllLayers && layer != 0))
                hx = exchangeVec[site - 1] * mxLHS + exchangeVec[site] * mxRHS + dipoleTerm + demagTerm;
        } else
            // All spins along x which are not within the driving region
            hx = exchangeVec[site - 1] * mxLHS + exchangeVec[site] * mxRHS + dipoleTerm + demagTerm;
    } else if (!isFm) {
        if (site >= drivingRegionLhs && site <= drivingRegionRhs) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (hasStaticDrive)
                hx = -1.0 * (exchangeVec[site - 1] * mxLHS + exchangeVec[site] * mxRHS + dynamicBiasField);
            else if (!hasStaticDrive)
                hx = -1.0 * (exchangeVec[site - 1] * mxLHS + exchangeVec[site] * mxRHS) + dynamicBiasField * cos(drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx = -1.0 * (exchangeVec[site - 1] * mxLHS + exchangeVec[site] * mxRHS);
    }

    return hx;
}
double EffectiveField::EffectiveFieldY(const int& site, const int& layer, const double& myLHS, const double& myMID, const double& myRHS,
                                                const double &dipoleTerm,  const double& demagTerm) {
    // The effective field (H_eff) y-component acting upon a given magnetic moment (site), abbreviated to 'hy'
    double hy;

    if (isFm) {
        hy = exchangeVec[site-1] * myLHS + exchangeVec[site] * myRHS + dipoleTerm + demagTerm;
    } else if (!isFm) {
        hy = -1.0 * (exchangeVec[site-1] * myLHS + exchangeVec[site] * myRHS);
    }

    return hy;
}
double EffectiveField::EffectiveFieldZ(const int& site, const int& layer, const double& mzLHS, const double& mzMID, const double& mzRHS,
                                                const double& dipoleTerm, const double& demagTerm) {
    // The effective field (H_eff) z-component acting upon a given magnetic moment (site), abbreviated to 'hz'
    double hz;

    if (isFm) {
        hz = exchangeVec[site-1] * mzLHS + exchangeVec[site] * mzRHS + dipoleTerm + demagTerm
                + GV.GetStaticBiasField();
    } else if (!isFm) {
        if (mzMID > 0)
            hz = GV.GetStaticBiasField() + anisotropyField - (exchangeVec[site-1] * mzLHS + exchangeVec[site] * mzRHS);
        else if (mzMID < 0)
            hz = GV.GetStaticBiasField() - anisotropyField - (exchangeVec[site-1] * mzLHS + exchangeVec[site] * mzRHS);
    }

    return hz;
}