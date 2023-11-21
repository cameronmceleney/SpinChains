//
// Created by Cameron McEleney on 31/10/2023.
//

#include "../include/EffectiveField.h"

EffectiveField::EffectiveField(SimulationParameters* sharedSimParams, 
                               SimulationStates* sharedSimStates, 
                               SimulationFlags* sharedSimFlags)
                                   
   : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}


double EffectiveField::EffectiveFieldX(const int& site, const int& layer, const double& mxLHS, const double& mxMID,
                                                const double& mxRHS, const double& dipoleTerm,
                                                const double& demagTerm, const double& current_time) {
    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hx;

    if (_simFlags->isFm) {
        if (site >= _simParams->drivingRegionLhs && site <= _simParams->drivingRegionRhs) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (_simFlags->driveAllLayers || layer == 0)
                hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm
                        + _simParams->dynamicBiasField * cos(_simParams->drivingAngFreq * current_time);
            else if  (_simFlags->hasStaticDrive)
                hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm
                        +_simParams->dynamicBiasField;
            else if  ((!_simFlags->driveAllLayers && layer != 0))
                hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm;
        } else
            // All spins along x which are not within the driving region
            hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dipoleTerm + demagTerm;
    } else if (!_simFlags->isFm) {
        if (site >= _simParams->drivingRegionLhs && site <= _simParams->drivingRegionRhs) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (_simFlags->hasStaticDrive)
                hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + _simParams->dynamicBiasField);
            else if (!_simFlags->hasStaticDrive)
                hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS) + _simParams->dynamicBiasField * cos(_simParams->drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS);
    }

    return hx;
}
double EffectiveField::EffectiveFieldY(const int& site, const int& layer, const double& myLHS, const double& myMID, const double& myRHS,
                                                const double &dipoleTerm,  const double& demagTerm) {
    // The effective field (H_eff) y-component acting upon a given magnetic moment (site), abbreviated to 'hy'
    double hy;

    if (_simFlags->isFm) {
        hy = _simStates->exchangeVec[site-1] * myLHS + _simStates->exchangeVec[site] * myRHS + dipoleTerm + demagTerm;
    } else if (!_simFlags->isFm) {
        hy = -1.0 * (_simStates->exchangeVec[site-1] * myLHS + _simStates->exchangeVec[site] * myRHS);
    }

    return hy;
}
double EffectiveField::EffectiveFieldZ(const int& site, const int& layer, const double& mzLHS, const double& mzMID, const double& mzRHS,
                                                const double& dipoleTerm, const double& demagTerm) {
    // The effective field (H_eff) z-component acting upon a given magnetic moment (site), abbreviated to 'hz'
    double hz;

    if (_simFlags->isFm) {
        hz = _simStates->exchangeVec[site-1] * mzLHS + _simStates->exchangeVec[site] * mzRHS + dipoleTerm + demagTerm
                + GV.GetStaticBiasField();
    } else if (!_simFlags->isFm) {
        if (mzMID > 0)
            hz = GV.GetStaticBiasField() + _simParams->anisotropyField - (_simStates->exchangeVec[site-1] * mzLHS + _simStates->exchangeVec[site] * mzRHS);
        else if (mzMID < 0)
            hz = GV.GetStaticBiasField() - _simParams->anisotropyField - (_simStates->exchangeVec[site-1] * mzLHS + _simStates->exchangeVec[site] * mzRHS);
    }

    return hz;
}