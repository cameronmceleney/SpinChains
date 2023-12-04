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

double
EffectiveFields::EffectiveFieldXClassic( const int &site, const int &layer, const double &mxLHS, const double &mxMID,
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
                hx += dipoleTerm + demagTerm +
                      _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * current_time);
            } else if ( _simFlags->isOscillatingZeemanStatic ) {
                hx = _simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS + dmiTerm;
                hx += +dipoleTerm + demagTerm + _simParams->oscillatingZeemanStrength;
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
            if ( _simFlags->isOscillatingZeemanStatic )
                hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS +
                             _simParams->oscillatingZeemanStrength);
            else if ( !_simFlags->isOscillatingZeemanStatic )
                hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS) +
                     _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx = -1.0 * (_simStates->exchangeVec[site - 1] * mxLHS + _simStates->exchangeVec[site] * mxRHS);
    }

    return hx;
}

double
EffectiveFields::EffectiveFieldYClassic( const int &site, const int &layer, const double &myLHS, const double &myMID,
                                         const double &myRHS, const double &dipoleTerm, const double &demagTerm,
                                         const double &dmiTerm ) {
    // The effective field (H_eff) y-component acting upon a given magnetic moment (site), abbreviated to 'hy'
    double hy;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        hy = _simStates->exchangeVec[site - 1] * myLHS + _simStates->exchangeVec[site] * myRHS - dmiTerm;
        hy += dipoleTerm + demagTerm;
    } else if ( !_simFlags->isFerromagnetic ) {
        hy = -1.0 * (_simStates->exchangeVec[site - 1] * myLHS + _simStates->exchangeVec[site] * myRHS);
    }

    return hy;
}

double
EffectiveFields::EffectiveFieldZClassic( const int &site, const int &layer, const double &mzLHS, const double &mzMID,
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

void
EffectiveFields::EffectiveFieldXTest( const int &currentSite, const int &layer, const std::vector<double> &mxTermsIn,
                                      std::vector<double> &effectiveFieldXTermsOut, const double &current_time ) {
    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hx = 0.0;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        if ( isSiteDriven(currentSite)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->shouldDriveAllLayers || layer == 0 ) {
                hx += _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
                hx += _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * current_time);
            } else if ( _simFlags->isOscillatingZeemanStatic ) {
                hx += _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
                hx += _simParams->oscillatingZeemanStrength;
            } else if ((!_simFlags->shouldDriveAllLayers && layer != 0)) {
                hx = _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                     + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
                hx += 0.0;
            }
        } else {
            // All spins along x which are not within the driving region
            hx += _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                  + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
            hx += 0.0;
        }
    } else if ( !_simFlags->isFerromagnetic ) {
        if ( isSiteDriven(currentSite)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->isOscillatingZeemanStatic )
                hx += -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1] +
                              _simParams->oscillatingZeemanStrength);
            else if ( !_simFlags->isOscillatingZeemanStatic )
                hx += -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1]) +
                      _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx += -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                          + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1]);
    }

    effectiveFieldXTermsOut[currentSite] = hx;
}

void
EffectiveFields::EffectiveFieldYTest( const int &currentSite, const int &layer, const std::vector<double> &myTermsIn,
                                      std::vector<double> &effectiveFieldYTermsOut ) {
    // The effective field (H_eff) y-component acting upon a given magnetic moment (site), abbreviated to 'hy'
    double hy = 0.0;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        hy += _simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
              + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1];
        hy += 0.0;
    } else if ( !_simFlags->isFerromagnetic ) {
        hy += -1.0 * (_simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1]);
    }

    effectiveFieldYTermsOut[currentSite] = hy;
}

void
EffectiveFields::EffectiveFieldZTest( const int &currentSite, const int &layer, const std::vector<double> &mzTermsIn,
                                      std::vector<double> &effectiveFieldZTermsOut ) {
    // The effective field (H_eff) z-component acting upon a given magnetic moment (site), abbreviated to 'hz'
    double hz = 0.0;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        hz += _simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
              + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1];
        hz += GV.GetStaticBiasField();
    } else if ( !_simFlags->isFerromagnetic ) {
        if ( mzTermsIn[currentSite] > 0 )
            hz += GV.GetStaticBiasField() + _simParams->anisotropyField -
                  (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                   + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
        else if ( mzTermsIn[currentSite] < 0 )
            hz += GV.GetStaticBiasField() - _simParams->anisotropyField -
                  (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                   + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
    }

    effectiveFieldZTermsOut[currentSite] = hz;
}

void EffectiveFields::EffectiveFieldsCombinedTest( const int &currentSite, const int &layer,
                                                   const std::vector<double> &mxTermsIn,
                                                   const std::vector<double> &myTermsIn,
                                                   const std::vector<double> &mzTermsIn,
                                                   std::vector<double> &effectiveFieldXTermsOut,
                                                   std::vector<double> &effectiveFieldYTermsOut,
                                                   std::vector<double> &effectiveFieldZTermsOut,
                                                   const double &currentTime ) {

    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    exchangeOut hExLocal;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        // hx terms
        if ( isSiteDriven(currentSite)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->shouldDriveAllLayers || layer == 0 ) {
                hExLocal.x += _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
                hExLocal.x += _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
            } else if ( _simFlags->isOscillatingZeemanStatic ) {
                hExLocal.x += _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
                hExLocal.x += _simParams->oscillatingZeemanStrength;
            } else if ((!_simFlags->shouldDriveAllLayers && layer != 0)) {
                hExLocal.x = _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                             + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
            }
        } else {
            // All spins along x which are not within the driving region
            hExLocal.x += _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                          + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1];
        }

        // hy terms
        hExLocal.y += _simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1];

        // hz terms
        hExLocal.z += _simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1];
        hExLocal.z += GV.GetStaticBiasField();

    } else if ( !_simFlags->isFerromagnetic ) {
        // hx terms
        if ( isSiteDriven(currentSite)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->isOscillatingZeemanStatic )
                hExLocal.x += -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                                      + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1] +
                                      _simParams->oscillatingZeemanStrength);
            else if ( !_simFlags->isOscillatingZeemanStatic )
                hExLocal.x += -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                                      + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1]) +
                              _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
        } else {
            // All spins along x which are not within the driving region
            hExLocal.x += -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                                  + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1]);
        }

        // hy terms
        hExLocal.y += -1.0 * (_simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1]);

        // hz terms
        if ( mzTermsIn[currentSite] > 0 )
            hExLocal.z += GV.GetStaticBiasField() + _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
        else if ( mzTermsIn[currentSite] < 0 )
            hExLocal.z += GV.GetStaticBiasField() - _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
    }

    if ( _simFlags->isFerromagnetic ) {
        hExLocal.y += _simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1];
    } else if ( !_simFlags->isFerromagnetic ) {
        hExLocal.y += -1.0 * (_simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1]);
    }

    if ( _simFlags->isFerromagnetic ) {
        hExLocal.z += _simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1];
        hExLocal.z += GV.GetStaticBiasField();
    } else if ( !_simFlags->isFerromagnetic ) {
        if ( mzTermsIn[currentSite] > 0 )
            hExLocal.z += GV.GetStaticBiasField() + _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
        else if ( mzTermsIn[currentSite] < 0 )
            hExLocal.z += GV.GetStaticBiasField() - _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
    }

    effectiveFieldXTermsOut[currentSite] = hExLocal.x;
    effectiveFieldYTermsOut[currentSite] = hExLocal.y;
    effectiveFieldZTermsOut[currentSite] = hExLocal.z;
}

std::array<double, 3> EffectiveFields::EffectiveFieldsCombinedTestExOnly( const int &currentSite, const int &layer,
                                                                          const std::vector<double> &mxTermsIn,
                                                                          const std::vector<double> &myTermsIn,
                                                                          const std::vector<double> &mzTermsIn ) {

    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hxLocal = 0.0, hyLocal = 0.0, hzLocal = 0.0;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        // hx terms
        hxLocal = _simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                    + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1]
                    -1.0 * _simParams->dmiConstant * (myTermsIn[currentSite + 1] - myTermsIn[currentSite - 1]);

        // hy terms
        hyLocal = _simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1]
                      + _simParams->dmiConstant * (mxTermsIn[currentSite + 1] - mxTermsIn[currentSite - 1]);

        // hz terms
        hzLocal = _simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]
                      + 0.0;

    } else if ( !_simFlags->isFerromagnetic ) {
        // hx terms
        hxLocal = -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTermsIn[currentSite - 1]
                    + _simStates->exchangeVec[currentSite] * mxTermsIn[currentSite + 1]);


        // hy terms
        hyLocal = -1.0 * (_simStates->exchangeVec[currentSite - 1] * myTermsIn[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * myTermsIn[currentSite + 1]);

        // hz terms
        if ( mzTermsIn[currentSite] > 0 )
            hzLocal = _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
        else if ( mzTermsIn[currentSite] < 0 )
            hzLocal =  -1.0 * _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTermsIn[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTermsIn[currentSite + 1]);
    }

    return {hxLocal, hyLocal, hzLocal};
}

std::array<double, 3> EffectiveFields::EffectiveFieldsCombinedTestDriveOnly( const int &currentSite, const int &layer,
                                                                             const std::vector<double> &mxTermsIn,
                                                                             const std::vector<double> &myTermsIn,
                                                                             const std::vector<double> &mzTermsIn,
                                                                             const double &currentTime ) {

    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hxLocal = 0.0, hyLocal = 0.0, hzLocal = 0.0;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        // hx terms
        if ( isSiteDriven(currentSite)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->shouldDriveAllLayers || layer == 0 ) {
                hxLocal = _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
            } else if ( _simFlags->isOscillatingZeemanStatic ) {
                hxLocal = _simParams->oscillatingZeemanStrength;
            }
        }

        // hy terms
        hyLocal = 0.0;

        // hz terms
        hzLocal = GV.GetStaticBiasField();

    } else if ( !_simFlags->isFerromagnetic ) {
        // hx terms
        if ( isSiteDriven(currentSite)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( _simFlags->isOscillatingZeemanStatic )
                hxLocal = _simParams->oscillatingZeemanStrength;
            else if ( !_simFlags->isOscillatingZeemanStatic )
                hxLocal = _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
        }

        // hy terms
        hyLocal = 0.0;

        // hz terms
        if ( mzTermsIn[currentSite] > 0 )
            hzLocal = GV.GetStaticBiasField();
        else if ( mzTermsIn[currentSite] < 0 )
            hzLocal = GV.GetStaticBiasField();
    }

    return {hxLocal, hyLocal, hzLocal};
}