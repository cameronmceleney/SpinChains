//
// Created by Cameron McEleney on 01/12/2023.
//

// Corresponding header
#include "../include/BiasFields.h"


BiasFields::BiasFields( SimulationParameters *sharedSimParams,
                        SimulationStates *sharedSimStates,
                        SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}

void BiasFields::calculateOneDimension( const int &currentLayer, const double &currentTime,
                                        std::vector<double> &biasFieldXOut, std::vector<double> &biasFieldYOut,
                                        std::vector<double> &biasFieldZOut ) {
    // This function is used for sequential calculations. Useful in small systems or when H_ext is complex
    std::array<double, 3> tempResultsContainer{};
    for ( int i = 1; i <= _simParams->systemTotalSpins; i++ ) {
        // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
        tempResultsContainer = _calculateBiasField1D(i, currentLayer, currentTime);
        biasFieldXOut[i] = tempResultsContainer[0];
        biasFieldYOut[i] = tempResultsContainer[1];
        biasFieldZOut[i] = tempResultsContainer[2];
    }
}

void BiasFields::calculateOneDimension( const int &currentLayer, const double &currentTime,
                                        std::vector<double> &biasFieldXOut, std::vector<double> &biasFieldYOut,
                                        std::vector<double> &biasFieldZOut, const bool &shouldUseTBB ) {
    // This function is used for parallel calculations. Useful in large systems or when H_ext is complex
    if ( shouldUseTBB ) {
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
                [&]( const tbb::blocked_range<int> tbbRange ) {
                    for ( int site = tbbRange.begin(); site <= tbbRange.end(); site++ ) {
                        // Need local vector to hold results to ensure this is threadsafe
                        std::array<double, 3> tempResultsLocal{0.0, 0.0, 0.0};
                        //tempResultsLocal = _calculateBiasField1D(i, currentLayer, currentTime, shouldUseTBB);
                        if ( site >= _simParams->drivingRegionLhs && site <= _simParams->drivingRegionRhs ) {
                            biasFieldXOut[site] += _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
                        } else {
                            biasFieldXOut[site] += 0.0;
                        }
                        biasFieldYOut[site] += 0.0;
                        biasFieldZOut[site] += GV.GetStaticBiasField();
                    }
        }, tbb::auto_partitioner());
    } else {
        throw std::invalid_argument("calculateOneDimension for exchange fields hasn't got CUDA implementation yet");
    }
}

std::array<double, 3>
BiasFields::_calculateBiasField1D( const int &currentSite, const int &currentLayer, const double &currentTime ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This version is verbose for easy debugging such as during the use of 'classic' methods
    // Only works when considering nearest-neighbours (NN)

    std::array<double, 3> biasFieldTerms{0.0, 0.0, 0.0};
    if ( _simFlags->isFerromagnetic ) {
        if ( _hasOscillatingZeeman(currentSite)) {
            // Currently assumes that the oscillating field will only be applied along x-axis
            if ( _simFlags->shouldDriveAllLayers || currentLayer == 0 )
                biasFieldTerms[0] +=
                        _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
            else if ( _simFlags->isOscillatingZeemanStatic )
                biasFieldTerms[0] += _simParams->oscillatingZeemanStrength;
        }
        if ( _simFlags->hasStaticZeeman && _simFlags->preferredDirection == 2 ) {
            // Currently assumes that the static field is uniform across z-axis
            biasFieldTerms[_simFlags->preferredDirection] += _simParams->staticZeemanStrength;
        }
    } else if ( !_simFlags->isFerromagnetic ) {
        if ( _hasOscillatingZeeman(currentSite)) {
            // Currently assumes that the oscillating field will only be applied along x-axis
            if ( _simFlags->shouldDriveAllLayers || currentLayer == 0 )
                biasFieldTerms[0] +=
                        _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
            else if ( _simFlags->isOscillatingZeemanStatic )
                biasFieldTerms[0] += _simParams->oscillatingZeemanStrength;
        }
        if ( _simFlags->hasStaticZeeman && _simFlags->preferredDirection == 2 ) {
            // Currently assumes that the static field is uniform across z-axis
            biasFieldTerms[_simFlags->preferredDirection] +=
                    _simParams->staticZeemanStrength;
        }
    }

    return biasFieldTerms;
}

std::array<double, 3>
BiasFields::_calculateBiasField1D( const int &currentSite, const int &currentLayer, const double &currentTime,
                                   const bool &shouldUseTBB ) {

    // TODO. This is a temp polymorphic version of _calculateDMIField1D that is threadsafe. Needs refinement
    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    /*
     * Code is currently hard-coded to only use nearest-neighbours (NN) heisenberg exchange interactions. But the
     * method is laid out this way to allow easier extension to include higher-order NN exchange interactions in
     * the future
     */

    if ( shouldUseTBB ) {
        // Haven't found a more efficient way to do this yet
        std::array<double, 3> biasFieldTerms{0.0, 0.0, 0.0};

        if ( _simFlags->isFerromagnetic ) {
            if ( _hasOscillatingZeeman(currentSite)) {
                // Currently assumes that the oscillating field will only be applied along x-axis
                if ( _simFlags->shouldDriveAllLayers || currentLayer == 0 )
                    biasFieldTerms[0] += _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
                else if ( _simFlags->isOscillatingZeemanStatic )
                    biasFieldTerms[0] += _simParams->oscillatingZeemanStrength;
            }
            if ( _simFlags->hasStaticZeeman ) {
                // Currently assumes that the static field is uniform across z-axis
                biasFieldTerms[2] += _simParams->staticZeemanStrength;
            }
        } else if ( !_simFlags->isFerromagnetic ) {
            if ( _hasOscillatingZeeman(currentSite)) {
                // Currently assumes that the oscillating field will only be applied along x-axis
                if ( _simFlags->shouldDriveAllLayers || currentLayer == 0 )
                    biasFieldTerms[0] +=
                            _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
                else if ( _simFlags->isOscillatingZeemanStatic )
                    biasFieldTerms[0] += _simParams->oscillatingZeemanStrength;
            }
            if ( _simFlags->hasStaticZeeman && _simFlags->preferredDirection == 2 ) {
                // Currently assumes that the static field is uniform across z-axis
                biasFieldTerms[_simFlags->preferredDirection] +=
                        _simParams->staticZeemanStrength;
            }
        }

        return biasFieldTerms;
    } else {
        throw std::invalid_argument("_calculateBiasFields1D hasn't got CUDA implementation yet");
    }
}

bool BiasFields::_hasOscillatingZeeman( const int &site ) {
    if ( _simFlags->shouldDriveDiscreteSites ) {
        for ( const int &discreteSite: _simStates->discreteDrivenSites )
            if ( site == discreteSite ) { return true; }
    }

    if ( site >= _simParams->drivingRegionLhs && site <= _simParams->drivingRegionRhs ) { return true; }

    // If no condition is met, then the site is not driven
    return false;
}

