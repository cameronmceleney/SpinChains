//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/EffectiveFields.h"

EffectiveFields::EffectiveFields( SimulationParameters *sharedSimParams,
                                  SimulationStates *sharedSimStates,
                                  SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}

bool EffectiveFields::isSiteDriven( const int &site ) {
    if ( _simFlags->shouldDriveDiscreteSites ) {
        for ( const int &discreteSite: _simStates->discreteDrivenSites )
            if ( site == discreteSite ) { return true; }
    }

    if ( site >= _simParams->drivingRegionLhs && site <= _simParams->drivingRegionRhs ) { return true; }

    // If no condition is met, then the site is not driven
    return false;
}

double EffectiveFields::EffectiveFieldX( const int &site, const int &layer, const double &mxLHS, const double &mxMID,
                                         const double &mxRHS, const double &dipoleTerm, const double &demagTerm,
                                         const double &dmiTerm, const double &current_time ) {
    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hx;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        if ( isSiteDriven(site)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->shouldDriveAllLayers || layer == 0 ) {
                hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dmiTerm;
                hx += dipoleTerm + demagTerm + _simParams->dynamicBiasField * cos(_simParams->drivingAngFreq * current_time);
            } else if ( _simFlags->isDriveStatic ) {
                hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dmiTerm;
                hx += + dipoleTerm + demagTerm + _simParams->dynamicBiasField;
            } else if ((!_simFlags->shouldDriveAllLayers && layer != 0)) {
                hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dmiTerm;
                hx += dipoleTerm + demagTerm;
            }
        } else {
            // All spins along x which are not within the driving region
            hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dmiTerm;
            hx += dipoleTerm + demagTerm;
        }
    } else if ( !_simFlags->isFerromagnetic ) {
        if ( isSiteDriven(site)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->isDriveStatic )
                hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS +
                             _simParams->dynamicBiasField);
            else if ( !_simFlags->isDriveStatic )
                hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS) +
                     _simParams->dynamicBiasField * cos(_simParams->drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS);
    }

    return hx;
}

double EffectiveFields::EffectiveFieldY( const int &site, const int &layer, const double &myLHS, const double &myMID,
                                         const double &myRHS, const double &dipoleTerm, const double &demagTerm,
                                         const double &dmiTerm ) {
    // The effective field (H_eff) y-component acting upon a given magnetic moment (site), abbreviated to 'hy'
    double hy;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        hy = _simStates->exchangeVec[site - 1] * myLHS + _simStates->exchangeVec[site] * myRHS - dmiTerm;
        hy += dipoleTerm + demagTerm;
    }
    else if ( !_simFlags->isFerromagnetic ) {
        hy = -1.0 * (_simStates->exchangeVec[site - 1] * myLHS + _simStates->exchangeVec[site] * myRHS);
    }

    return hy;
}

double EffectiveFields::EffectiveFieldZ( const int &site, const int &layer, const double &mzLHS, const double &mzMID,
                                         const double &mzRHS,
                                         const double &dipoleTerm, const double &demagTerm ) {
    // The effective field (H_eff) z-component acting upon a given magnetic moment (site), abbreviated to 'hz'
    double hz;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic )
        hz = _simStates->exchangeVec[site - 1] * mzLHS + _simStates->exchangeVec[site] * mzRHS + dipoleTerm +
             demagTerm + GV.GetStaticBiasField();
    else if ( !_simFlags->isFerromagnetic ) {
        if ( mzMID > 0 )
            hz = GV.GetStaticBiasField() + _simParams->anisotropyField -
                 (_simStates->exchangeVec[site - 1] * mzLHS + _simStates->exchangeVec[site] * mzRHS);
        else if ( mzMID < 0 )
            hz = GV.GetStaticBiasField() - _simParams->anisotropyField -
                 (_simStates->exchangeVec[site - 1] * mzLHS + _simStates->exchangeVec[site] * mzRHS);
    }

    return hz;
}