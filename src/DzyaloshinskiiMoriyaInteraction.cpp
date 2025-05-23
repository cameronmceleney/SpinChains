//
// Created by Cameron McEleney on 22/11/2023.
//

// Corresponding header
#include "../include/DzyaloshinskiiMoriyaInteraction.h"

DzyaloshinskiiMoriyaInteraction::DzyaloshinskiiMoriyaInteraction( SimulationParameters *sharedSimParams,
                                                                  SimulationStates *sharedSimStates,
                                                                  SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {

    // Currently assuming that the DMI vector is along the z-axis (i.e. choosing 'z' over 'x' when symmetry breaking is
    // along the 'y' axis)

}

CommonStructures::Vector3D
DzyaloshinskiiMoriyaInteraction::calculateClassic( const int &currentSite, const CommonStructures::Vector2D &mxTerms,
                                                   const CommonStructures::Vector2D &myTerms,
                                                   const CommonStructures::Vector2D &mzTerms ) {
    // Only use for debugging!!
    return _calculateDMIFieldClassic(currentSite, mxTerms, myTerms, mzTerms);
}

void DzyaloshinskiiMoriyaInteraction::calculateOneDimension( const std::vector<double> &mxTerms,
                                                             const std::vector<double> &myTerms,
                                                             const std::vector<double> &mzTerms,
                                                             std::vector<double> &dmiXOut,
                                                             std::vector<double> &dmiYOut,
                                                             std::vector<double> &dmiZOut ) {

    for ( int i = 1; i <= _simParams->systemTotalSpins; i++ ) {
        // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
        auto tempResultsContainer = _calculateDMIField1D(i, mxTerms, myTerms, mzTerms);
        dmiXOut[i] = tempResultsContainer.x();
        dmiYOut[i] = tempResultsContainer.y();
        dmiZOut[i] = tempResultsContainer.z();
    }
}

void DzyaloshinskiiMoriyaInteraction::calculateOneDimension( const std::vector<double> &mxTerms,
                                                             const std::vector<double> &myTerms,
                                                             const std::vector<double> &mzTerms,
                                                             std::vector<double> &dmiXOut,
                                                             std::vector<double> &dmiYOut, std::vector<double> &dmiZOut,
                                                             const bool &shouldUseTBB ) {
    // This function is used for parallel calculations. Useful in large systems or when DMI field is complex
    if ( shouldUseTBB ) {
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
                          [&]( const tbb::blocked_range<int> tbbRange ) {
                              for ( int i = tbbRange.begin(); i <= tbbRange.end(); i++ ) {
                                  // Need local vector to hold results to ensure this is threadsafe
                                  CommonStructures::Vector3D tempResultsLocal = _calculateDMIField1D(i, mxTerms, myTerms,
                                                                                                mzTerms, shouldUseTBB);
                                  dmiXOut[i] = tempResultsLocal.x();
                                  dmiYOut[i] = tempResultsLocal.y();
                                  dmiZOut[i] = tempResultsLocal.z();
                              }
                          });
    }
    else {
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
        auto tempResultsContainer = _calculateDMIField3D(i, mxTerms, myTerms, mzTerms);
        dmiXOut[i] = tempResultsContainer.x();
        dmiYOut[i] = tempResultsContainer.y();
        dmiZOut[i] = tempResultsContainer.z();
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
                                  CommonStructures::Vector3D tempResultsContainer = _calculateDMIField3D(i, mxTerms, myTerms,
                                                                                                    mzTerms,
                                                                                                    shouldUseTBB);
                                  dmiXOut[i] = tempResultsContainer.x();
                                  dmiYOut[i] = tempResultsContainer.y();
                                  dmiZOut[i] = tempResultsContainer.z();
                              }
                          });
    } else {
        throw std::invalid_argument("calculateThreeDimensions hasn't got CUDA implementation yet");
    }
}


CommonStructures::Vector3D
DzyaloshinskiiMoriyaInteraction::_calculateDMIField1D( const int &currentSite, const std::vector<double> &mxTerms,
                                                       const std::vector<double> &myTerms,
                                                       const std::vector<double> &mzTerms ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This is 1D so there is no DMI vector, only one non-zero component which is assumed to be along the z-axis

    CommonStructures::Vector3D tempResults(0.0, 0.0, 0.0);

    if ( currentSite <= _simParams->numSpinsInABC or
         currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
        // Guard clause to ensure that the current site doesn't lie within the aborbing boundary condition region
        return tempResults;
    }

    // As cross produt of vectors A x B = - B x A, we can use the same function for both cases. So only process leftwise
    CommonStructures::Vector3D originSite = {mxTerms[currentSite], myTerms[currentSite], mzTerms[currentSite]};
    CommonStructures::Vector3D influencingSite = {mxTerms[currentSite - 1], myTerms[currentSite - 1],
                                             mzTerms[currentSite - 1]};

    CommonStructures::Vector3D originCrossInfluencingSites = _crossProduct(originSite, influencingSite);
    // In 1D the typical _crossProduct(_dmiVector, originCrossInfluencingSites) simply becomes a dot product with only the z-component being non-zero

    for (int i = 0; i < 3; i++) {tempResults[i] = _simParams->dmiVector[i] * originCrossInfluencingSites[i];}
    return tempResults;
}

CommonStructures::Vector3D
DzyaloshinskiiMoriyaInteraction::_calculateDMIField1D( const int &currentSite, const std::vector<double> &mxTerms,
                                                       const std::vector<double> &myTerms,
                                                       const std::vector<double> &mzTerms, const bool &shouldUseTBB ) {

    // TODO. This is a temp polymorphic version of _calculateDMIField1D that is threadsafe. Needs refinement
    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This is 1D so there is no DMI vector, only one non-zero component which is assumed to be along the z-axis

    CommonStructures::Vector3D tempResults(0.0, 0.0, 0.0);

    if ( shouldUseTBB ) {
        // As cross product of vectors A x B = - B x A, we can use the same function for both cases. So only process leftwise
        // TODO. Find a way to remove as many vector creations as possible
        CommonStructures::Vector3D originSiteLocal = {mxTerms[currentSite], myTerms[currentSite], mzTerms[currentSite]};
        CommonStructures::Vector3D influencingSiteLocal = {mxTerms[currentSite - 1], myTerms[currentSite - 1],
                                                      mzTerms[currentSite - 1]};

        CommonStructures::Vector3D originCrossInfluencingSitesLocal = _crossProduct(originSiteLocal, influencingSiteLocal);

        // In 1D the typical _crossProduct(_dmiVector, originCrossInfluencingSites) simply becomes a dot product with only the z-component being non-zero
        // dmiDotSites = {0.0, 0.0, _simParams->dmiVector.z() * originCrossInfluencingSites[2]};
        for (int i = 0; i < 3; i++) {tempResults[i] = _simParams->dmiVector[i] * originCrossInfluencingSitesLocal[i];}

        return tempResults;
    }
    else { throw std::invalid_argument("_calculateDMIField1D hasn't got CUDA implementation yet"); }
}

CommonStructures::Vector3D
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
    CommonStructures::Vector3D originSite = {mxTerms[currentSite], myTerms[currentSite], mzTerms[currentSite]};
    CommonStructures::Vector3D influencingSite = {mxTerms[currentSite - 1], myTerms[currentSite - 1],
                                             mzTerms[currentSite - 1]};

    CommonStructures::Vector3D originCrossInfluencingSites = _crossProduct(originSite, influencingSite);
    return _crossProduct(_simParams->dmiVector, originCrossInfluencingSites);  // This is 'dmiCrossSites'
}

CommonStructures::Vector3D
DzyaloshinskiiMoriyaInteraction::_calculateDMIField3D( auto &currentSite, const std::vector<double> &mxTerms,
                                                       const std::vector<double> &myTerms,
                                                       const std::vector<double> &mzTerms,
                                                       const bool &shouldUseTBB ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // Method is written to be easily extended to 3D which will happen soon
    if ( shouldUseTBB ) {
        if ( currentSite <= _simParams->numSpinsInABC ||
             currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
            // Guard clause to ensure that the current site doesn't lie within the absorbing boundary condition region
            return {0.0, 0.0, 0.0};
        }

        // As cross product of vectors A x B = - B x A, we can use the same function for both cases. So only process leftwise
        CommonStructures::Vector3D originSite = {mxTerms[currentSite], myTerms[currentSite], mzTerms[currentSite]};
        CommonStructures::Vector3D influencingSite = {mxTerms[currentSite - 1], myTerms[currentSite - 1],
                                                 mzTerms[currentSite - 1]};

        CommonStructures::Vector3D originCrossInfluencingSites = _crossProduct(originSite, influencingSite);
        return _crossProduct(_simParams->dmiVector, originCrossInfluencingSites);  // This is 'dmiCrossSites'
    } else {
        throw std::invalid_argument("_calculateDMIField3D hasn't got CUDA implementation yet");
    }
}

CommonStructures::Vector3D
DzyaloshinskiiMoriyaInteraction::_calculateDMIFieldClassic( auto &currentSite, const CommonStructures::Vector2D &mxTerms,
                                                            const CommonStructures::Vector2D &myTerms,
                                                            const CommonStructures::Vector2D &mzTerms ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This is 1D so there is no DMI vector, only one non-zero component which is assumed to be along the z-axis

    // As cross product of vectors A x B = - B x A, we can use the same function for both cases. So only process leftwise
    CommonStructures::Vector3D originSite = {mxTerms[1], myTerms[1], mzTerms[1]};
    CommonStructures::Vector3D influencingSite = {mxTerms[0], myTerms[0], mzTerms[0]};

    CommonStructures::Vector3D originCrossInfluencingSites = _crossProduct(originSite, influencingSite);
    // In 1D the typical _crossProduct(_dmiVector, originCrossInfluencingSites) simply becomes a dot product with only the z-component being non-zero

    return {_simParams->dmiVector.x() * originCrossInfluencingSites[0],
            _simParams->dmiVector.y() * originCrossInfluencingSites[1],
            _simParams->dmiVector.z() * originCrossInfluencingSites[2]};  // This is 'dmiDotSites'
}

CommonStructures::Vector3D DzyaloshinskiiMoriyaInteraction::_crossProduct( const CommonStructures::Vector3D &iSite,
                                                                      const CommonStructures::Vector3D &jSite ) {

    CommonStructures::Vector3D crossProductVector{0.0, 0.0, 0.0};

    crossProductVector.x() =  iSite.y() * jSite.z() - iSite.z() * jSite.y();
    crossProductVector.y() = -iSite.x() * jSite.z() + iSite.z() * jSite.x();
    crossProductVector.z() =  iSite.x() * jSite.y() - iSite.y() * jSite.x();

    return crossProductVector;
}