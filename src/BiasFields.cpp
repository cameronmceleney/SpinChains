//
// Created by Cameron McEleney on 01/12/2023.
//

// Corresponding header
#include "../include/BiasFields.h"


BiasFields::BiasFields( SimulationParameters *sharedSimParams,
                        SimulationStates *sharedSimStates,
                        SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}

void BiasFields::calculateOneDimension( const int &currentLayer, const double &currentTime, const std::vector<double> &mzTermsIn,
                                        std::vector<double> &biasFieldXOut, std::vector<double> &biasFieldYOut,
                                        std::vector<double> &biasFieldZOut ) {
    // This function is used for sequential calculations. Useful in small systems or when H_ext is complex
    std::array<double, 3> tempResultsContainer{};
    bool shouldUseTBB = true;

    for ( int i = 1; i <= _simParams->systemTotalSpins; i++ ) {
        // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
        tempResultsContainer = _calculateBiasField1D(i, currentLayer, currentTime, mzTermsIn[i], shouldUseTBB);
        biasFieldXOut[i] = tempResultsContainer[0];
        biasFieldYOut[i] = tempResultsContainer[1];
        biasFieldZOut[i] = tempResultsContainer[2];
    }
}

void BiasFields::calculateOneDimension( const int &currentLayer, const double &currentTime,
                                        const std::vector<double> &mzTermsIn, std::vector<std::atomic<double>> &biasFieldXOut,
                                        std::vector<std::atomic<double>> &biasFieldYOut, std::vector<std::atomic<double>> &biasFieldZOut,
                                        const bool &shouldUseTBB ) {
    // This function is used for parallel calculations. Useful in large systems or when H_ext is complex

    if ( shouldUseTBB ) {

        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
            [&](const tbb::blocked_range<int>& range) {
                for (int site = range.begin(); site < range.end(); site++) {
                    std::array<double, 3> tempBiasResults = _calculateBiasField1D( site, currentLayer, currentTime, mzTermsIn[site], shouldUseTBB);

                    biasFieldXOut[site].fetch_add(tempBiasResults[0]);
                    biasFieldYOut[site].fetch_add(tempBiasResults[1]);
                    biasFieldZOut[site].fetch_add(tempBiasResults[2]);
                }
        }, tbb::auto_partitioner());

    } else {
        throw std::invalid_argument("calculateOneDimension for exchange fields hasn't got CUDA implementation yet");
    }
}

void BiasFields::calculateOneDimension( const int &currentLayer, const double &currentTime,
                                        const std::vector<double> &mzTermsIn, std::vector<double> &biasFieldXOut,
                                        std::vector<double> &biasFieldYOut, std::vector<double> &biasFieldZOut,
                                        const bool &shouldUseTBB ) {
    // Doesn't work right now

    if ( shouldUseTBB ) {
        /*
        // Thread-local storage for each component
        tbb::combinable<std::vector<double>> localFieldX([&]{ return std::vector<double>(_simParams->systemTotalSpins + 2, 0.0); });
        tbb::combinable<std::vector<double>> localFieldY([&]{ return std::vector<double>(_simParams->systemTotalSpins + 2, 0.0); });
        tbb::combinable<std::vector<double>> localFieldZ([&]{ return std::vector<double>(_simParams->systemTotalSpins + 2, 0.0); });

        // Parallel computation
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
            [&](const tbb::blocked_range<int>& range) {
                auto& localX = localFieldX.local();
                auto& localY = localFieldY.local();
                auto& localZ = localFieldZ.local();

                for (int i = range.begin(); i <= range.end(); i++) {
                    std::array<double, 3> fieldContribution = _calculateBiasField1D( i, currentLayer, currentTime, mzTermsIn[i], shouldUseTBB);
                    localX[i] += fieldContribution[0];
                    localY[i] += fieldContribution[1];
                    localZ[i] += fieldContribution[2];
                }
            }, tbb::auto_partitioner());

        // Combining the results directly into the output fields
        localFieldX.combine_each([&](const std::vector<double>& v) {
            for (size_t i = 1; i <= _simParams->systemTotalSpins; i++) biasFieldXOut[i] += v[i];
        });
        localFieldY.combine_each([&](const std::vector<double>& v) {
            for (size_t i = 1; i <= _simParams->systemTotalSpins; i++) biasFieldYOut[i] += v[i];
        });
        localFieldZ.combine_each([&](const std::vector<double>& v) {
            for (size_t i = 1; i <= _simParams->systemTotalSpins; i++) biasFieldZOut[i] += v[i];
        });
        */
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
            [&](const tbb::blocked_range<int>& range) {
                for (int site = range.begin(); site <= range.end(); site++) {
                    std::array<double, 3> tempBiasResults = _calculateBiasField1D( site, currentLayer, currentTime, mzTermsIn[site], shouldUseTBB);

                    std::atomic_ref<double>(biasFieldXOut[site]).fetch_add(tempBiasResults[0]);
                    std::atomic_ref<double>(biasFieldYOut[site]).fetch_add(tempBiasResults[1]);
                    std::atomic_ref<double>(biasFieldZOut[site]).fetch_add(tempBiasResults[2]);
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
                biasFieldTerms[0] =
                        _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
            else if ( _simFlags->isOscillatingZeemanStatic )
                biasFieldTerms[0] = _simParams->oscillatingZeemanStrength;
        }
        if ( _simFlags->hasStaticZeeman && _simFlags->preferredDirection == 2 ) {
            // Currently assumes that the static field is uniform across z-axis
            biasFieldTerms[_simFlags->preferredDirection] = _simParams->staticZeemanStrength;
        }
    } else if ( !_simFlags->isFerromagnetic ) {
        if ( _hasOscillatingZeeman(currentSite)) {
            // Currently assumes that the oscillating field will only be applied along x-axis
            if ( _simFlags->shouldDriveAllLayers || currentLayer == 0 )
                biasFieldTerms[0] =
                        _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
            else if ( _simFlags->isOscillatingZeemanStatic )
                biasFieldTerms[0] = _simParams->oscillatingZeemanStrength;
        }
        if ( _simFlags->hasStaticZeeman && _simFlags->preferredDirection == 2 ) {
            // Currently assumes that the static field is uniform across z-axis
            biasFieldTerms[_simFlags->preferredDirection] =
                    _simParams->staticZeemanStrength;
        }
    }

    return biasFieldTerms;
}

std::array<double, 3>
BiasFields::_calculateBiasField1D( const int &currentSite, const int &currentLayer, const double &currentTime,
                                   const double &mzTermAtSite, const bool &shouldUseTBB ) {

    // TODO. This is a temp polymorphic version of _calculateDMIField1D that is threadsafe. Needs refinement
    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    /*
     * Code is currently hard-coded to only use nearest-neighbours (NN) heisenberg exchange interactions. But the
     * method is laid out this way to allow easier extension to include higher-order NN exchange interactions in
     * the future
     */

    if ( shouldUseTBB ) {
        // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
        double hxLocal = 0.0, hyLocal = 0.0, hzLocal = 0.0;

        // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

        if ( _simFlags->isFerromagnetic ) {
            // hx terms
            if ( _hasOscillatingZeeman(currentSite)) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                if ( _simFlags->shouldDriveAllLayers || currentLayer == 0 ) {
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
            if ( _hasOscillatingZeeman(currentSite)) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                if ( _simFlags->isOscillatingZeemanStatic )
                    hxLocal = _simParams->oscillatingZeemanStrength;
                else if ( !_simFlags->isOscillatingZeemanStatic )
                    hxLocal = _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
            }

            // hy terms
            hyLocal = 0.0;

            // hz terms
            if ( mzTermAtSite > 0 )
                hzLocal = GV.GetStaticBiasField();
            else if ( mzTermAtSite < 0 )
                hzLocal = GV.GetStaticBiasField();
        }

        return {hxLocal, hyLocal, hzLocal};
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

