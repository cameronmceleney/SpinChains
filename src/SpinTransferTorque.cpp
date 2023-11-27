//
// Created by Cameron McEleney on 27/11/2023.
//

// Corresponding header
#include "../include/SpinTransferTorque.h"

SpinTransferTorque::SpinTransferTorque( SimulationParameters *sharedSimParams,
                                        SimulationStates *sharedSimStates,
                                        SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {
}

void
SpinTransferTorque::calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                           const std::vector<double> &mzTerms, std::vector<double> &sttXOut,
                                           std::vector<double> &sttYOut, std::vector<double> &sttZOut ) {


    for ( int iSite = 1; iSite <= _simParams->systemTotalSpins; iSite++ ) {
        if ( iSite <= _simParams->numSpinsInABC || iSite > (_simParams->numSpinsInChain + _simParams->numSpinsInABC)) {
            // Check boundary conditions
            sttXOut[iSite] = 0.0;
            sttYOut[iSite] = 0.0;
            sttZOut[iSite] = 0.0;
        }

        std::array<double, 3> spCurrentDensity = _calculateExchangeDrivenSpinCurrentDensity(iSite, mxTerms, myTerms,
                                                                                            mzTerms);

        // Calculate the spin current. Might need to be more complex
        std::array<double, 3> spinCurrent = {_simParams->spinPolarisation * spCurrentDensity[0],
                                             _simParams->spinPolarisation * spCurrentDensity[1],
                                             _simParams->spinPolarisation * spCurrentDensity[2]};

        std::array<double, 3> currentMTerms = {mxTerms[iSite], myTerms[iSite], mzTerms[iSite]};

        // Calculate the adiabatic and non-adiabatic terms. Are these both unique to STT, or am I recalculating the heisenberg exchange/dipolar fields etc?
        std::array<double, 3> adiabaticTerm = _crossProduct(currentMTerms, spinCurrent);
        std::array<double, 3> nonAdiabaticTerm = _crossProduct(currentMTerms, adiabaticTerm);

        // Return STT elements so other methods can add them to the LLG equation. Unsure if I need both the adiabatic/non-adiabatic?
        sttXOut[iSite] = _simParams->gyroMagConst * adiabaticTerm[0] +
                         _simParams->spinTransferEfficiency * _simParams->gilbertDamping *
                         nonAdiabaticTerm[0];
        sttYOut[iSite] = _simParams->gyroMagConst * adiabaticTerm[1] +
                         _simParams->spinTransferEfficiency * _simParams->gilbertDamping *
                         nonAdiabaticTerm[1];
        sttZOut[iSite] = _simParams->gyroMagConst * adiabaticTerm[2] +
                         _simParams->spinTransferEfficiency * _simParams->gilbertDamping *
                         nonAdiabaticTerm[2];
    }
}

std::array<double, 3> SpinTransferTorque::_calculateExchangeDrivenSpinCurrentDensity( const int &currentSite,
                                                                                      const std::vector<double> &mxTerms,
                                                                                      const std::vector<double> &myTerms,
                                                                                      const std::vector<double> &mzTerms ) {
    // Initialize the spin current density vector
    std::array<double, 3> spin_current_density = {0.0, 0.0, 0.0};
    if ( currentSite <= _simParams->numSpinsInABC ||
         currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
        // Check boundary conditions
        return {0.0, 0.0, 0.0};
    }

    if ( _simStates->exchangeVec[currentSite] == 0 ) {
        // Guard clause to ensure that the exchange vector is not zero
        return {0.0, 0.0, 0.0};
    }

    double latticeConstant = std::sqrt(_simParams->exchangeStiffness / _simStates->exchangeVec[currentSite]);

    if ( std::isinf(latticeConstant)) {
        // Guard clause to ensure that the lattice constant is not infinite (backup test / temporary)
        throw std::runtime_error(std::string("Lattice constant is infinite!"));
    }

    // All checks passed; proceed with finding exchange driven spin current density
    std::array<double, 3> originSite = {mxTerms[currentSite], myTerms[currentSite], mzTerms[currentSite]};
    std::array<double, 3> influencingSite = {mxTerms[currentSite - 1], myTerms[currentSite - 1],
                                             mzTerms[currentSite - 1]};

    // Calculate the gradient of magnetization (central difference)
    std::array<double, 3> grad_m;
    for ( int i = 0; i < 3; ++i ) {
        grad_m[i] = (originSite[i] - influencingSite[i]) / (2 * latticeConstant);
    }

    // Calculate the exchange-driven spin current density
    for ( int i = 0; i < 3; ++i ) {
        spin_current_density[i] = -_simParams->exchangeStiffness * grad_m[i];
    }

    // For the boundaries, one may use forward or backward difference
    // or apply periodic boundary conditions depending on the physical system.

    return spin_current_density;
}

std::array<double, 3> SpinTransferTorque::_crossProduct( const std::array<double, 3> &iSite,
                                                         const std::array<double, 3> &jSite ) {
    if ( iSite.size() != 3 || jSite.size() != 3 )
        throw std::invalid_argument("One or more input vectors to DMI::crossProduct() are not size 3");

    std::array<double, 3> crossProductVector{};
    crossProductVector[0] = iSite[1] * jSite[2] - iSite[2] * jSite[1];
    crossProductVector[1] = -iSite[0] * jSite[2] + iSite[2] * jSite[0];
    crossProductVector[2] = iSite[0] * jSite[1] - iSite[1] * jSite[0];

    return crossProductVector;
}

std::array<double, 3> SpinTransferTorque::_crossProduct( const std::array<double, 3> &iSite,
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