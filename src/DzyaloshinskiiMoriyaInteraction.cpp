//
// Created by Cameron McEleney on 22/11/2023.
//

// Corresponding header
#include "../include/DzyaloshinskiiMoriyaInteraction.h"

DzyaloshinskiiMoriyaInteraction::DzyaloshinskiiMoriyaInteraction( SimulationParameters *sharedSimParams,
                                                                  SimulationStates *sharedSimStates,
                                                                  SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {

    _dmiVector = {0, 0, _simParams->dmiConstant};
}

std::vector<double>
DzyaloshinskiiMoriyaInteraction::calculateClassic( const int &currentSite, const std::vector<double> &mxTerms,
                                                   const std::vector<double> &myTerms,
                                                   const std::vector<double> &mzTerms ) {
    // Only use for debugging!!
    return _calculateDMIFieldClassic(currentSite, mxTerms, myTerms, mzTerms);
}

void DzyaloshinskiiMoriyaInteraction::calculateOneDimension( const std::vector<double> &mxTerms,
                                                             const std::vector<double> &myTerms,
                                                             const std::vector<double> &mzTerms,
                                                             std::vector<double> &dmiXOut,
                                                             std::vector<double> &dmiYOut,
                                                             std::vector<double> &dmiZOut ) {
    // This function is used for sequential calculations. Useful in small systems or when DMI field is simple
    for ( int i = 1; i <= _simParams->systemTotalSpins; i++ ) {
        // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
        _tempResultsContainer = _calculateDMIField1D(i, mxTerms, myTerms, mzTerms);
        dmiXOut[i] = _tempResultsContainer[0];
        dmiYOut[i] = _tempResultsContainer[1];
        dmiZOut[i] = _tempResultsContainer[2];
    }
}

void DzyaloshinskiiMoriyaInteraction::calculateOneDimension( const std::vector<double> &mxTerms,
                                                             const std::vector<double> &myTerms,
                                                             const std::vector<double> &mzTerms,
                                                             std::vector<double> &dmiXOut,
                                                             std::vector<double> &dmiYOut, std::vector<double> &dmiZOut,
                                                             bool shouldUseTBB ) {
    // This function is used for parallel calculations. Useful in large systems or when DMI field is complex
    if ( shouldUseTBB ) {
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
                          [&]( const tbb::blocked_range<int> tbbRange ) {
                              for ( int i = tbbRange.begin(); i <= tbbRange.end(); i++ ) {
                                  // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
                                  _tempResultsContainer = _calculateDMIField1D(i, mxTerms, myTerms, mzTerms);
                                  dmiXOut[i] = _tempResultsContainer[0];
                                  dmiYOut[i] = _tempResultsContainer[1];
                                  dmiZOut[i] = _tempResultsContainer[2];
                              }
                          });
    } else {
        throw std::invalid_argument("calculateOneDimension hasn't got CUDA implementation yet");
    }
}

void DzyaloshinskiiMoriyaInteraction::calculateThreeDimensions( const std::vector<double> &mxTerms,
                                                                const std::vector<double> &myTerms,
                                                                const std::vector<double> &mzTerms,
                                                                std::vector<double> &dmiXOut,
                                                                std::vector<double> &dmiYOut,
                                                                std::vector<double> &dmiZOut ) {
    // This function is used for sequential calculations. Useful in small systems or when DMI field is simple
    for ( int i = 1; i <= _simParams->systemTotalSpins; i++ ) {
        // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
        _tempResultsContainer = _calculateDMIField3D(i, mxTerms, myTerms, mzTerms);
        dmiXOut[i] = _tempResultsContainer[0];
        dmiYOut[i] = _tempResultsContainer[1];
        dmiZOut[i] = _tempResultsContainer[2];
    }
}

void DzyaloshinskiiMoriyaInteraction::calculateThreeDimensions( const std::vector<double> &mxTerms,
                                                                const std::vector<double> &myTerms,
                                                                const std::vector<double> &mzTerms,
                                                                std::vector<double> &dmiXOut,
                                                                std::vector<double> &dmiYOut,
                                                                std::vector<double> &dmiZOut,
                                                                bool shouldUseTBB ) {
    // This function is used for parallel calculations. Useful in large systems or when DMI field is complex

    if ( shouldUseTBB ) {
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
                          [&]( const tbb::blocked_range<int> tbbRange ) {
                              for ( int i = tbbRange.begin(); i <= tbbRange.end(); i++ ) {
                                  // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
                                  _tempResultsContainer = _calculateDMIField3D(i, mxTerms, myTerms, mzTerms);
                                  dmiXOut[i] = _tempResultsContainer[0];
                                  dmiYOut[i] = _tempResultsContainer[1];
                                  dmiZOut[i] = _tempResultsContainer[2];
                              }
                          });
    } else {
        throw std::invalid_argument("calculateOneDimension hasn't got CUDA implementation yet");
    }
}


std::vector<double>
DzyaloshinskiiMoriyaInteraction::_calculateDMIField1D( auto &currentSite, const std::vector<double> &mxTerms,
                                                       const std::vector<double> &myTerms,
                                                       const std::vector<double> &mzTerms ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This is 1D so there is no DMI vector, only one non-zero component which is assumed to be along the z-axis

    if ( currentSite <= _simParams->numSpinsInABC ||
         currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
        // Guard clause to ensure that the current site doesn't lie within the aborbing boundary condition region
        return {0.0, 0.0, 0.0};
    }

    // As cross produt of vectors A x B = - B x A, we can use the same function for both cases. So only process leftwise
    // TODO. Find a way to remove as many vector creations as possible
    _originSite = {mxTerms[currentSite], myTerms[currentSite], mzTerms[currentSite]};
    _influencingSite = {mxTerms[currentSite - 1], myTerms[currentSite - 1], mzTerms[currentSite - 1]};

    _originCrossInfluencingSites = _crossProduct(_originSite, _influencingSite);
    // In 1D the typical _crossProduct(_dmiVector, originCrossInfluencingSites) simply becomes a dot product with only the z-component being non-zero
    // dmiDotSites = {0.0, 0.0, _simParams->dmiConstant * originCrossInfluencingSites[2]};
    return {0.0, 0.0, _simParams->dmiConstant * _originCrossInfluencingSites[2]};
}

std::vector<double>
DzyaloshinskiiMoriyaInteraction::_calculateDMIField3D( auto &currentSite, const std::vector<double> &mxTerms,
                                                       const std::vector<double> &myTerms,
                                                       const std::vector<double> &mzTerms ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // Method is written to be easily extended to 3D which will happen soon

    if ( currentSite <= _simParams->numSpinsInABC ||
         currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
        // Guard clause to ensure that the current site doesn't lie within the aborbing boundary condition region
        return {0.0, 0.0, 0.0};
    }

    // As cross produt of vectors A x B = - B x A, we can use the same function for both cases. So only process leftwise
    _originSite = {mxTerms[currentSite], myTerms[currentSite], mzTerms[currentSite]};
    _influencingSite = {mxTerms[currentSite - 1], myTerms[currentSite - 1],
                                           mzTerms[currentSite - 1]};

    _originCrossInfluencingSites = _crossProduct(_originSite, _influencingSite);
    return _crossProduct(_dmiVector, _originCrossInfluencingSites);  // This is 'dmiCrossSites'
}

std::vector<double>
DzyaloshinskiiMoriyaInteraction::_calculateDMIFieldClassic( auto &currentSite, const std::vector<double> &mxTerms,
                                                            const std::vector<double> &myTerms,
                                                            const std::vector<double> &mzTerms ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This is 1D so there is no DMI vector, only one non-zero component which is assumed to be along the z-axis

    if ( currentSite <= _simParams->numSpinsInABC ||
         currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
        // Guard clause to ensure that the current site doesn't lie within the aborbing boundary condition region
        return {0.0, 0.0, 0.0};
    }

    // As cross produt of vectors A x B = - B x A, we can use the same function for both cases. So only process leftwise
    _originSite = {mxTerms[1], myTerms[1], mzTerms[1]};
    _influencingSite = {mxTerms[0], myTerms[0], mzTerms[0]};

    _originCrossInfluencingSites = _crossProduct(_originSite, _influencingSite);
    // In 1D the typical _crossProduct(_dmiVector, originCrossInfluencingSites) simply becomes a dot product with only the z-component being non-zero
    return {0.0, 0.0, _simParams->dmiConstant * _originCrossInfluencingSites[2]};  // This is 'dmiDotSites'
}

std::vector<double> DzyaloshinskiiMoriyaInteraction::_crossProduct( const std::vector<double> &iSite,
                                                                    const std::vector<double> &jSite ) {
    // Method is currently unoptimised for debugging purposes. Will use an inline method in the future
    // TODO. This will likely be a performance bottleneck in the future
    if ( iSite.size() != 3 || jSite.size() != 3 )
        throw std::invalid_argument("One or more input vectors to DMI::crossProduct() are not size 3");

    std::vector<double> crossProductVector(3);
    crossProductVector[0] = iSite[1] * jSite[2] - iSite[2] * jSite[1];
    crossProductVector[1] = -iSite[0] * jSite[2] + iSite[2] * jSite[0];
    crossProductVector[2] = iSite[0] * jSite[1] - iSite[1] * jSite[0];

    return crossProductVector;
}