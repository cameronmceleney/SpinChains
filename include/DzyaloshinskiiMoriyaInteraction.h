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
                                std::vector<double> &dmiYOut, std::vector<double> &dmiZOut, const bool &shouldUseTBB );

    void calculateThreeDimensions( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                   const std::vector<double> &mzTerms, std::vector<double> &dmiXOut,
                                   std::vector<double> &dmiYOut, std::vector<double> &dmiZOut );

    void calculateThreeDimensions( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                   const std::vector<double> &mzTerms, std::vector<double> &dmiXOut,
                                   std::vector<double> &dmiYOut, std::vector<double> &dmiZOut, bool shouldUseTBB );

    std::array<double, 3> calculateClassic( const int &currentSite, const std::array<double, 2> &mxTerms,
                                            const std::array<double, 2> &myTerms,
                                            const std::array<double, 2> &mzTerms );

private:
    SimulationParameters *_simParams; // Non-owning pointer to SimulationParameters
    SimulationStates *_simStates;
    SimulationFlags *_simFlags;
private:
    // Empty contains to be constants reused throughout the component's lifetime instead recreating new each method call
    std::array<double, 3> _dmiVector{};

    /**
     * To be used in all cases where readability is key for debugging, and all single-threaded cases
     * @param iSite
     * @param jSite
     * @return
     */
    inline std::array<double, 3> _crossProduct( const std::array<double, 3> &iSite, const std::array<double, 3> &jSite );

    /**
     * To be used in multi-threaded cases and is optimised for efficiency
     * @param iSite
     * @param jSite
     * @return
     */
    inline std::array<double, 3> _crossProduct( const std::array<double, 3> &iSite, const std::array<double, 3> &jSite,
                                                const bool & shouldUseTBB);

    /**
     * To be used in all cases where readability is key for debugging, and all single-threaded cases
     * @param iSite
     * @param jSite
     * @return
     */
    inline std::array<double, 3> _dotProduct( const std::array<double, 3> &iSite, const std::array<double, 3> &jSite );

    /**
     * To be used in multi-threaded cases and is optimised for efficiency
     * @param iSite
     * @param jSite
     * @return
     */
    inline std::array<double, 3> _dotProduct( const std::array<double, 3> &iSite, const std::array<double, 3> &jSite,
                                                const bool & shouldUseTBB);

    std::array<double, 3> _calculateDMIField1D( const int &currentSite, const std::vector<double> &mxTerms,
                                              const std::vector<double> &myTerms, const std::vector<double> &mzTerms );

    std::array<double, 3>
    _calculateDMIField1D( const int &currentSite, const std::vector<double> &mxTerms,
                          const std::vector<double> &myTerms,
                          const std::vector<double> &mzTerms, const bool &shouldUseTBB );

    std::array<double, 3> _calculateDMIField3D( auto &currentSite, const std::vector<double> &mxTerms,
                                              const std::vector<double> &myTerms, const std::vector<double> &mzTerms );
    std::array<double, 3> _calculateDMIField3D( auto &currentSite, const std::vector<double> &mxTerms,
                                          const std::vector<double> &myTerms, const std::vector<double> &mzTerms,
                                          const bool &shouldUseTBB);

    std::array<double, 3>
    _calculateDMIFieldClassic( auto &currentSite, const std::array<double, 2> &mxTerms,
                               const std::array<double, 2> &myTerms,
                               const std::array<double, 2> &mzTerms );
};


#endif //SPINCHAINS_DZYALOSHINSKIIMORIYAINTERACTION_H
