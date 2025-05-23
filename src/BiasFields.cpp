//
// Created by Cameron McEleney on 01/12/2023.
//

// Corresponding header
#include "../include/BiasFields.h"


BiasFields::BiasFields( SimulationParameters *sharedSimParams,
                        SimulationStates *sharedSimStates,
                        SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}


// Refactored code
// 1D public
void BiasFields::calculateOneDimension(const int &currentLayer, const double &currentTime,
                                       const std::vector<double> &mzTermsIn, std::vector<double> &biasFieldXOut,
                                       std::vector<double> &biasFieldYOut, std::vector<double> &biasFieldZOut,
                                       const CommonStructures::Parallelisations &parallelisationMode) {
    calculateBiasField(1, currentLayer, currentTime, mzTermsIn, biasFieldXOut, biasFieldYOut, biasFieldZOut, false, parallelisationMode);
}

void BiasFields::calculateOneDimension(const int &currentLayer, const double &currentTime,
                                       const std::vector<double> &mzTermsIn, std::vector<std::atomic<double>> &biasFieldXOut,
                                       std::vector<std::atomic<double>> &biasFieldYOut,
                                       std::vector<std::atomic<double>> &biasFieldZOut,
                                       const CommonStructures::Parallelisations &parallelisationMode) {
    calculateBiasField(1, currentLayer, currentTime, mzTermsIn, biasFieldXOut, biasFieldYOut, biasFieldZOut, true, parallelisationMode);
}

template <typename T>
void BiasFields::calculateBiasField(int dimension, const int &currentLayer, const double &currentTime,
                                    const std::vector<double> &mzTermsIn, std::vector<T> &biasFieldXOut,
                                    std::vector<T> &biasFieldYOut, std::vector<T> &biasFieldZOut,
                                    bool isAtomic, CommonStructures::Parallelisations parallelFlag) {
    if (parallelFlag == CommonStructures::Parallelisations::Multithreaded) {
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins + 1),
            [&](const tbb::blocked_range<int> &range) {
                for (int site = range.begin(); site < range.end(); site++) {
                    auto tempBiasResults = calculateBiasField(dimension, site, currentLayer, currentTime, mzTermsIn);
                    addBiasField(site, tempBiasResults, biasFieldXOut, biasFieldYOut, biasFieldZOut, isAtomic);
                }
            }, tbb::auto_partitioner());
    }
    else if (parallelFlag == CommonStructures::Parallelisations::Sequential) {
        for (int site = 1; site <= _simParams->systemTotalSpins; site++) {
            auto tempBiasResults = calculateBiasField(dimension, site, currentLayer, currentTime, mzTermsIn);
            addBiasField(site, tempBiasResults, biasFieldXOut, biasFieldYOut, biasFieldZOut, isAtomic);
        }
    }
    else {
        throw std::invalid_argument("BiasFields::calculateBiasField: Invalid parallelisation flag");
    }
}

CommonStructures::Vector3D BiasFields::calculateBiasField(int dimension, int site, const int &currentLayer, const double &currentTime,
                                                     const std::vector<double> &mzTermsIn) {

    CommonStructures::Vector3D biasFieldTerms{0.0, 0.0, 0.0};

    // No need for scalingFactor as a variable while we're only making a single check
    //double scalingFactor = 1.0; // Implies no scaling changes are made to the bias field

    if (dimension == 1) { biasFieldTerms = _calculateBiasField1D(site, currentLayer, currentTime); }
    else { throw std::invalid_argument("BiasFields::calculateBiasField: Dimensionality beyond 1 not supported"); }

    // If this flag is not enabled, we have an early return opportunity
    if (not _simFlags->hasGradientRegionForOscillatingZeeman) { return biasFieldTerms; }

    // Check driving region map to see if current site has a scaling factor
    auto it = _simStates->dRGradientMap.find(site);
    if ( it != _simStates->dRGradientMap.end()) {
        // Only change the x-component of the bias field as this is the only component that will be driven dynamically
        biasFieldTerms.x() *= it->second;
    }

    return biasFieldTerms;
}

CommonStructures::Vector3D
BiasFields::_calculateBiasField1D( int site, const int &currentLayer, const double &currentTime ) {
    CommonStructures::Vector3D biasFieldTerms{0.0, 0.0, 0.0};

    if (_hasOscillatingZeeman(site)) {
        // Currently assumes that the oscillating field will only be applied along x-axis
        if (_simFlags->shouldDriveAllLayers || currentLayer == 0) {
            biasFieldTerms[Axis::x] += _simParams->oscillatingZeemanStrength * cos(_simParams->drivingAngFreq * currentTime);
        }
        else if (_simFlags->isOscillatingZeemanStatic) {
            biasFieldTerms[Axis::x] += _simParams->oscillatingZeemanStrength;
        }
    }

    if (_simFlags->hasStaticZeeman) {
        biasFieldTerms[_simFlags->equilibriumOrientation] += _simParams->staticZeemanStrength[_simFlags->equilibriumOrientation];
    }

    /*
     * Could add an IF statement to check `simFlags->_simFlags->isFerromagnetic` in the future
     */

    return biasFieldTerms;
}

template <typename T>
void BiasFields::addBiasField(int site, const CommonStructures::Vector3D &biasField, std::vector<T> &biasFieldXOut,
                              std::vector<T> &biasFieldYOut, std::vector<T> &biasFieldZOut, bool isAtomic) {
    biasFieldXOut[site] += biasField[Axis::x];
    biasFieldYOut[site] += biasField[Axis::y];
    biasFieldZOut[site] += biasField[Axis::z];
}

// Specialization for std::atomic<double>
template <>
void BiasFields::addBiasField<std::atomic<double>>(int site, const CommonStructures::Vector3D &biasField, std::vector<std::atomic<double>> &biasFieldXOut,
                                                   std::vector<std::atomic<double>> &biasFieldYOut, std::vector<std::atomic<double>> &biasFieldZOut, bool isAtomic) {
    biasFieldXOut[site].fetch_add(biasField[Axis::x], std::memory_order_relaxed);
    biasFieldYOut[site].fetch_add(biasField[Axis::y], std::memory_order_relaxed);
    biasFieldZOut[site].fetch_add(biasField[Axis::z], std::memory_order_relaxed);
    // std::atomic_ref<double>(biasFieldXOut[site]).fetch_add(tempBiasResults[0]);
    // std::atomic_ref<double>(biasFieldYOut[site]).fetch_add(tempBiasResults[1]);
    // std::atomic_ref<double>(biasFieldZOut[site]).fetch_add(tempBiasResults[2]);
}

bool BiasFields::_hasOscillatingZeeman( const int &site ) {

    if (_simFlags->shouldDriveDiscreteSites) {
        for (const int &discreteSite: _simStates->discreteDrivenSites) {
            if (site == discreteSite) { return true; }
        }
    }

    // Check simple bounds second; possible for user to set these alongside `_simFlags->shouldDriveDiscreteSites`
    // which leads to logic errors
    if (site >= _simParams->drivingRegion.lhsSite && site <= _simParams->drivingRegion.rhsSite) { return true; }

    return false;
}

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
            CommonStructures::Vector3D fieldContribution = _calculateBiasField1D( i, currentLayer, currentTime, mzTermsIn[i], shouldUseTBB);
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
