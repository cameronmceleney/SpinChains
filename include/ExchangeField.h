//
// Created by Cameron McEleney on 01/12/2023.
//

#ifndef SPINCHAINS_EXCHANGEFIELD_H
#define SPINCHAINS_EXCHANGEFIELD_H

// C++ Standard Library
#include <vector>

// C++ Third Party Libraries
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/mutex.h>

// C++ User Libraries (General)
#include "CommonLibs.h"
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"

class ExchangeField {
private:
    SimulationParameters *_simParams; // Non-owning pointer to SimulationParameters
    SimulationStates *_simStates;
    SimulationFlags *_simFlags;

private:

    std::array<double, 3> _calculateExchangeField1D( const int &currentSite, const std::vector<double> &mxTerms,
                                                     const std::vector<double> &myTerms,
                                                     const std::vector<double> &mzTerms );

    std::array<double, 3>
    _calculateExchangeField1D( const int &currentSite, const std::vector<double> &mxTerms,
                               const std::vector<double> &myTerms,
                               const std::vector<double> &mzTerms, const bool &shouldUseTBB );

        // Empty contains to be constants reused throughout the component's lifetime instead recreating new each method call
    std::array<double, 3> _dmiVector{};

    /**
     * To be used in all cases where readability is key for debugging, and all single-threaded cases
     * @param iSite
     * @param jSite
     * @return
     */
    static inline std::array<double, 3> _crossProduct( const std::array<double, 3> &iSite, const std::array<double, 3> &jSite );

    /**
     * To be used in multi-threaded cases and is optimised for efficiency
     * @param iSite
     * @param jSite
     * @return
     */
    static inline std::array<double, 3> _crossProduct( const std::array<double, 3> &iSite, const std::array<double, 3> &jSite,
                                                const bool & shouldUseTBB);


    std::array<double, 3> _calculateDMI1D( const int &currentSite, const std::vector<double> &mxTerms,
                                              const std::vector<double> &myTerms, const std::vector<double> &mzTerms );

    std::array<double, 3>
    _calculateDMI1D( const int &currentSite, const std::vector<double> &mxTerms,
                          const std::vector<double> &myTerms,
                          const std::vector<double> &mzTerms, const bool &shouldUseTBB );

public:
    explicit ExchangeField( SimulationParameters *sharedSimParams,
                            SimulationStates *sharedSimStates,
                            SimulationFlags *sharedSimFlags );

    ~ExchangeField() = default;

public:
    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<double> &exchangeXOut,
                                std::vector<double> &exchangeYOut, std::vector<double> &exchangeZOut );

    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<std::atomic<double>> &exchangeXOut,
                                std::vector<std::atomic<double>> &exchangeYOut, std::vector<std::atomic<double>> &exchangeZOut,
                                const bool &shouldUseTBB );

    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<double> &exchangeXOut,
                                std::vector<double> &exchangeYOut, std::vector<double> &exchangeZOut,
                                const bool &shouldUseTBB );
};


#endif //SPINCHAINS_EXCHANGEFIELD_H
