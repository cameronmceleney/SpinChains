//
// Created by Cameron McEleney on 27/11/2023.
//

#ifndef SPINCHAINS_SPINTRANSFERTORQUE_H
#define SPINCHAINS_SPINTRANSFERTORQUE_H


// C++ Standard Library
#include <vector>

// C++ Third Party Libraries
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

// C++ User Libraries (General)
#include "../libs/CommonDefinitions.h"
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"

class SpinTransferTorque {
public:
    explicit SpinTransferTorque( SimulationParameters *sharedSimParams,
                                              SimulationStates *sharedSimStates,
                                              SimulationFlags *sharedSimFlags );

    ~SpinTransferTorque() = default;

public:
    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<double> &sttXOut,
                                std::vector<double> &sttYOut, std::vector<double> &sttZOut );

private:
    SimulationParameters *_simParams; // Non-owning pointer to SimulationParameters
    SimulationStates *_simStates;
    SimulationFlags *_simFlags;
private:
    // Empty contains to be constants reused throughout the component's lifetime instead recreating new each method call
    CommonStructures::Vector3D _calculateExchangeDrivenSpinCurrentDensity( const int &currentSite, const std::vector<double> &mxTerms,
                                                                      const std::vector<double> &myTerms,
                                                                      const std::vector<double> &mzTerms);
    CommonStructures::Vector3D _calculateSTT1D();

    /**
     * To be used in all cases where readability is key for debugging, and all single-threaded cases
     * @param iSite
     * @param jSite
     * @return
     */
    inline CommonStructures::Vector3D _crossProduct( const CommonStructures::Vector3D &iSite, const CommonStructures::Vector3D &jSite );

    /**
     * To be used in multi-threaded cases and is optimised for efficiency
     * @param iSite
     * @param jSite
     * @return
     */
    inline CommonStructures::Vector3D _crossProduct( const CommonStructures::Vector3D &iSite, const CommonStructures::Vector3D &jSite,
                                                const bool & shouldUseTBB);
};


#endif //SPINCHAINS_SPINTRANSFERTORQUE_H
