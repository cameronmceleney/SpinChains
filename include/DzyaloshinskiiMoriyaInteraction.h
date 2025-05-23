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
#include "../libs/CommonDefinitions.h"
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

    CommonStructures::Vector3D calculateClassic( const int &currentSite, const CommonStructures::Vector2D &mxTerms,
                                            const CommonStructures::Vector2D &myTerms,
                                            const CommonStructures::Vector2D &mzTerms );

private:
    SimulationParameters *_simParams; // Non-owning pointer to SimulationParameters
    SimulationStates *_simStates;
    SimulationFlags *_simFlags;
private:
    // Empty contains to be constants reused throughout the component's lifetime instead recreating new each method call
    /**
     * To be used in all cases where readability is key for debugging, and all single-threaded cases
     * @param iSite
     * @param jSite
     * @return
     */
    static inline CommonStructures::Vector3D _crossProduct( const CommonStructures::Vector3D &iSite, const CommonStructures::Vector3D &jSite );

    /**
     * To be used in all cases where readability is key for debugging, and all single-threaded cases
     * @param iSite
     * @param jSite
     * @return
     */
    inline CommonStructures::Vector3D _dotProduct( const CommonStructures::Vector3D &iSite, const CommonStructures::Vector3D &jSite );

    /**
     * To be used in multi-threaded cases and is optimised for efficiency
     * @param iSite
     * @param jSite
     * @return
     */
    inline CommonStructures::Vector3D _dotProduct( const CommonStructures::Vector3D &iSite, const CommonStructures::Vector3D &jSite,
                                                const bool & shouldUseTBB);

    CommonStructures::Vector3D _calculateDMIField1D( const int &currentSite, const std::vector<double> &mxTerms,
                                              const std::vector<double> &myTerms, const std::vector<double> &mzTerms );

    CommonStructures::Vector3D
    _calculateDMIField1D( const int &currentSite, const std::vector<double> &mxTerms,
                          const std::vector<double> &myTerms,
                          const std::vector<double> &mzTerms, const bool &shouldUseTBB );

    CommonStructures::Vector3D _calculateDMIField3D( auto &currentSite, const std::vector<double> &mxTerms,
                                              const std::vector<double> &myTerms, const std::vector<double> &mzTerms );
    CommonStructures::Vector3D _calculateDMIField3D( auto &currentSite, const std::vector<double> &mxTerms,
                                          const std::vector<double> &myTerms, const std::vector<double> &mzTerms,
                                          const bool &shouldUseTBB);

    CommonStructures::Vector3D
    _calculateDMIFieldClassic( auto &currentSite, const CommonStructures::Vector2D &mxTerms,
                               const CommonStructures::Vector2D &myTerms,
                               const CommonStructures::Vector2D &mzTerms );
};


#endif //SPINCHAINS_DZYALOSHINSKIIMORIYAINTERACTION_H
