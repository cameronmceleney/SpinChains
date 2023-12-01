//
// Created by Cameron McEleney on 01/12/2023.
//

// Corresponding header
#include "../include/ExchangeField.h"

ExchangeField::ExchangeField( SimulationParameters *sharedSimParams,
                              SimulationStates *sharedSimStates,
                              SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {

    _dmiVector = {0, 0, _simParams->dmiConstant};
}

void ExchangeField::calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                           const std::vector<double> &mzTerms, std::vector<double> &exchangeXOut,
                                           std::vector<double> &exchangeYOut, std::vector<double> &exchangeZOut ) {
    // This function is used for sequential calculations. Useful in small systems or when H_ex is simple
    std::array<double, 3> tempResultsContainerHeisenberg{};
    std::array<double, 3> tempResultsContainerDMI{};

    for ( int i = 1; i <= _simParams->systemTotalSpins; i++ ) {
        // Used boundary limits [1, _simParams->systemTotalSpins] (inclusive) is intentional
        tempResultsContainerHeisenberg = _calculateExchangeField1D(i, mxTerms, myTerms, mzTerms);
        tempResultsContainerDMI = _calculateDMI1D(i, mxTerms, myTerms, mzTerms);

        exchangeXOut[i] += tempResultsContainerHeisenberg[0] + tempResultsContainerDMI[0];
        exchangeYOut[i] += tempResultsContainerHeisenberg[1] + tempResultsContainerDMI[1];
        exchangeZOut[i] += tempResultsContainerHeisenberg[2] + tempResultsContainerDMI[2];
    }
}

void ExchangeField::calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                           const std::vector<double> &mzTerms, std::vector<double> &exchangeXOut,
                                           std::vector<double> &exchangeYOut, std::vector<double> &exchangeZOut,
                                           const bool &shouldUseTBB ) {
    // This function is used for parallel calculations. Useful in large systems or when H_ex is complex
    if ( shouldUseTBB ) {
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
                          [&]( const tbb::blocked_range<int> tbbRange ) {
                              for ( int i = tbbRange.begin(); i <= tbbRange.end(); i++ ) {
                                  // Need local vector to hold results to ensure this is threadsafe
                                  std::array<double, 3> tempResultsLocal = _calculateExchangeField1D(i, mxTerms,
                                                                                                     myTerms,
                                                                                                     mzTerms,
                                                                                                     shouldUseTBB);
                                  exchangeXOut[i] = tempResultsLocal[0];
                                  exchangeYOut[i] = tempResultsLocal[1];
                                  exchangeZOut[i] = tempResultsLocal[2];

                                  // Overwrite tempResultsLocal with DMI results to save memory
                                  tempResultsLocal = _calculateDMI1D(i, mxTerms, myTerms, mzTerms, shouldUseTBB);
                                  exchangeXOut[i] += tempResultsLocal[0];
                                  exchangeYOut[i] += tempResultsLocal[1];
                                  exchangeZOut[i] += tempResultsLocal[2];
                              }
                          });
    } else {
        throw std::invalid_argument("calculateOneDimension for exchange fields hasn't got CUDA implementation yet");
    }
}

std::array<double, 3>
ExchangeField::_calculateExchangeField1D( const int &currentSite, const std::vector<double> &mxTerms,
                                          const std::vector<double> &myTerms,
                                          const std::vector<double> &mzTerms ) {

    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This version is verbose for easy debugging such as during the use of 'classic' methods
    // Only works when considering nearest-neighbours (NN)

    // Separate defintions for LHS and RHS to minimise changes of a mistake occurring
    double exchangeLhs = _simStates->exchangeVec[currentSite - 1];
    double exchangeRhs = _simStates->exchangeVec[currentSite];

    std::array<double, 3> heisenbergExchangeTerms{0.0, 0.0, 0.0};
    if ( _simFlags->isFerromagnetic ) {
        heisenbergExchangeTerms[0] = exchangeLhs * mxTerms[currentSite - 1] +
                                     exchangeRhs * mxTerms[currentSite + 1];
        heisenbergExchangeTerms[1] = exchangeLhs * myTerms[currentSite - 1] +
                                     exchangeRhs * myTerms[currentSite + 1];

        heisenbergExchangeTerms[2] = exchangeLhs * mzTerms[currentSite - 1] +
                                     exchangeRhs * mzTerms[currentSite + 1];
    } else if ( !_simFlags->isFerromagnetic ) {
        heisenbergExchangeTerms[0] = -1.0 * (exchangeLhs * mxTerms[currentSite - 1] +
                                             exchangeRhs * mxTerms[currentSite + 1]);
        heisenbergExchangeTerms[1] = -1.0 * (exchangeLhs * myTerms[currentSite - 1] +
                                             exchangeRhs * myTerms[currentSite + 1]);

        if ( mzTerms[currentSite] > 0 )
            // TODO. Anisotropy field needs moved to its own class in the future. Written this way to be explicit
            heisenbergExchangeTerms[2] = +1.0 * _simParams->anisotropyField -
                                         (exchangeLhs * mzTerms[currentSite - 1] +
                                          exchangeRhs * mzTerms[currentSite + 1]);
        else if ( mzTerms[currentSite] < 0 )
            heisenbergExchangeTerms[2] = -1.0 * _simParams->anisotropyField -
                                         (exchangeLhs * mzTerms[currentSite - 1] +
                                          exchangeRhs * mzTerms[currentSite + 1]);
    }

    return heisenbergExchangeTerms;
}

std::array<double, 3>
ExchangeField::_calculateExchangeField1D( const int &currentSite, const std::vector<double> &mxTerms,
                                          const std::vector<double> &myTerms,
                                          const std::vector<double> &mzTerms, const bool &shouldUseTBB ) {

    // TODO. This is a temp polymorphic version of _calculateDMIField1D that is threadsafe. Needs refinement
    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    /*
     * Code is currently hard-coded to only use nearest-neighbours (NN) heisenberg exchange interactions. But the
     * method is laid out this way to allow easier extension to include higher-order NN exchange interactions in
     * the future
     */

    if ( shouldUseTBB && _simFlags->isFerromagnetic ) {
            return {
                    _simStates->exchangeVec[currentSite - 1] * mxTerms[currentSite - 1] +
                    _simStates->exchangeVec[currentSite] * mxTerms[currentSite + 1],
                    _simStates->exchangeVec[currentSite - 1] * myTerms[currentSite - 1] +
                    _simStates->exchangeVec[currentSite] * myTerms[currentSite + 1],
                    _simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1] +
                    _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1],
            };
    } else
        throw std::invalid_argument("_calculateExchangeField1D hasn't got CUDA implementation yet");
}

std::array<double, 3>
ExchangeField::_calculateDMI1D( const int &currentSite, const std::vector<double> &mxTerms,
                                const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms ) {

    /* Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
     * Methodically finds general solution to H_DMI = D_{i-1, i} \cdot (m_{i-1} \times m_{i}) + D_{i, i+1} \cdot (m_{i} \times m_{i+1}
     * for when there are only two nearest neighbours
     *
     * Uses Eq. 3 from https://doi.org/10.1103/PhysRevB.107.224421 to explicitly write the return statements for the case
     * where DMI only involves two nearest neighbours (NN).
     */

    std::array<double, 3> siteLhs = {mxTerms[currentSite - 1], myTerms[currentSite - 1], mzTerms[currentSite - 1]};
    std::array<double, 3> siteRhs = {mxTerms[currentSite + 1], myTerms[currentSite + 1],
                                     mzTerms[currentSite + 1]};

    std::array<double, 3> dmiEquation = {  _dmiVector[1] * (siteRhs[2] - siteLhs[2]) - _dmiVector[2] * (siteRhs[1] - siteLhs[1]),
                                          -_dmiVector[0] * (siteRhs[2] - siteLhs[2]) + _dmiVector[2] * (siteRhs[0] - siteLhs[0]),
                                           _dmiVector[0] * (siteRhs[1] - siteLhs[1]) - _dmiVector[1] * (siteRhs[0] - siteLhs[0])
                                        };

    // In 1D the typical _crossProduct(_dmiVector, originCrossInfluencingSites) simply becomes a dot product with only the z-component being non-zero
    return {dmiEquation[0], dmiEquation[1], dmiEquation[2]};
}

std::array<double, 3>
ExchangeField::_calculateDMI1D( const int &currentSite, const std::vector<double> &mxTerms,
                                const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, const bool &shouldUseTBB ) {

    /* TODO. This is a temp polymorphic version of _calculateDMIField1D that is threadsafe. Needs refinement
     * Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
     *
     * Quickly finds solution to H_DMI = D_{i-1, i} \cdot (m_{i-1} \times m_{i}) + D_{i, i+1} \cdot (m_{i} \times m_{i+1}
     * by exploiting how only D_z is non-zero in this 1D system
     *
     * Uses Eq. 3 from https://doi.org/10.1103/PhysRevB.107.224421 to explicitly write the return statements for the case
     * where DMI only involves two nearest neighbours (NN).
     */


    if ( shouldUseTBB ) {
        return {-1.0 * _simParams->dmiConstant * (myTerms[currentSite + 1] - myTerms[currentSite + 1]),
                _simParams->dmiConstant * (mxTerms[currentSite + 1] - mxTerms[currentSite + 1]),
                0.0};
    } else {
        throw std::invalid_argument("_calculateDMIField1D hasn't got CUDA implementation yet");
    }
}

std::array<double, 3> ExchangeField::_crossProduct( const std::array<double, 3> &iSite,
                                                    const std::array<double, 3> &jSite ) {
    if ( iSite.size() != 3 || jSite.size() != 3 )
        throw std::invalid_argument("One or more input vectors to DMI::crossProduct() are not size 3");

    std::array<double, 3> crossProductVector{};
    crossProductVector[0] = iSite[1] * jSite[2] - iSite[2] * jSite[1];
    crossProductVector[1] = -iSite[0] * jSite[2] + iSite[2] * jSite[0];
    crossProductVector[2] = iSite[0] * jSite[1] - iSite[1] * jSite[0];

    return crossProductVector;
}

std::array<double, 3> ExchangeField::_crossProduct( const std::array<double, 3> &iSite,
                                                    const std::array<double, 3> &jSite,
                                                    const bool &shouldUseTBB ) {
    // All needed tests are done by calling method

    // Able to return result directly; no need to use temp containers; wasted time and memory on heap
    if ( shouldUseTBB )
        return {
                iSite[1] * jSite[2] - iSite[2] * jSite[1],
                -iSite[0] * jSite[2] + iSite[2] * jSite[0],
                iSite[0] * jSite[1] - iSite[1] * jSite[0]
        };
    else
        throw std::invalid_argument("_crossProduct hasn't got CUDA implementation yet");
}
