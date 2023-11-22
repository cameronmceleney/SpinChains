//
// Created by Cameron McEleney on 22/11/2023.
//

#ifndef SPINCHAINS_DZYALOSHINSKIIMORIYAINTERACTION_H
#define SPINCHAINS_DZYALOSHINSKIIMORIYAINTERACTION_H

// C++ Standard Library
#include <vector>

// C++ Third Party Libraries
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

// C++ User Libraries (General)
#include "CommonLibs.h"
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"

class DzyaloshinskiiMoriyaInteraction {
public:
    explicit DzyaloshinskiiMoriyaInteraction( SimulationParameters *sharedSimParams,
                                              SimulationStates *sharedSimStates,
                                              SimulationFlags *sharedSimFlags );

    ~DzyaloshinskiiMoriyaInteraction() = default;

public:
    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<double> &dmiXOut,
                                std::vector<double> &dmiYOut, std::vector<double> &dmiZOut );

    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<double> &dmiXOut,
                                std::vector<double> &dmiYOut, std::vector<double> &dmiZOut, bool shouldUseTBB );

    void calculateThreeDimensions( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                   const std::vector<double> &mzTerms, std::vector<double> &dmiXOut,
                                   std::vector<double> &dmiYOut, std::vector<double> &dmiZOut );

    void calculateThreeDimensions( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                   const std::vector<double> &mzTerms, std::vector<double> &dmiXOut,
                                   std::vector<double> &dmiYOut, std::vector<double> &dmiZOut, bool shouldUseTBB );

    std::vector<double> calculateClassic( const int &currentSite, const std::vector<double> &mxTerms,
                                          const std::vector<double> &myTerms,
                                          const std::vector<double> &mzTerms );

private:
    SimulationParameters *_simParams; // Non-owning pointer to SimulationParameters
    SimulationStates *_simStates;
    SimulationFlags *_simFlags;
private:
    // Empty contains to be constants reused throughout the component's lifetime instead recreating new each method call
    std::vector<double> _tempResultsContainer = {0.0, 0.0, 0.0};
    std::vector<double> _originSite = {0.0, 0.0, 0.0};
    std::vector<double> _influencingSite = {0.0, 0.0, 0.0};
    std::vector<double> _originCrossInfluencingSites = {0.0, 0.0, 0.0};
    std::vector<double> _dmiVector;


    std::vector<double> _crossProduct( const std::vector<double> &iSite, const std::vector<double> &jSite );

    std::vector<double> _calculateDMIField1D( auto &currentSite, const std::vector<double> &mxTerms,
                                              const std::vector<double> &myTerms, const std::vector<double> &mzTerms );

    std::vector<double> _calculateDMIField3D( auto &currentSite, const std::vector<double> &mxTerms,
                                              const std::vector<double> &myTerms, const std::vector<double> &mzTerms );

    std::vector<double>
    _calculateDMIFieldClassic( auto &currentSite, const std::vector<double> &mxTerms,
                               const std::vector<double> &myTerms,
                               const std::vector<double> &mzTerms );
};


#endif //SPINCHAINS_DZYALOSHINSKIIMORIYAINTERACTION_H
