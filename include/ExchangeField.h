//
// Created by Cameron McEleney on 01/12/2023.
//

#ifndef SPINCHAINS_EXCHANGEFIELD_H
#define SPINCHAINS_EXCHANGEFIELD_H

// C++ Standard Library
#include <array>
#include <vector>

// C++ Third Party Libraries
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/mutex.h>

// C++ User Libraries (General)
#include "../libs/CommonDefinitions.h"
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
public:
    explicit ExchangeField( SimulationParameters *sharedSimParams,
                            SimulationStates *sharedSimStates,
                            SimulationFlags *sharedSimFlags );

    ~ExchangeField() = default;

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
     * To be used in multi-threaded cases and is optimised for efficiency
     * @param iSite
     * @param jSite
     * @return
     */

private:
    template <typename T>
    void calculateExchangeField( int dimension, const std::vector<double> &mxTerms,
                                 const std::vector<double> &myTerms, const std::vector<double> &mzTerms,
                                 std::vector<T> &exchangeXOut, std::vector<T> &exchangeYOut,
                                 std::vector<T> &exchangeZOut,
                                 CommonStructures::Parallelisations parallelFlag );

    CommonStructures::Vector3D _calculateHeisenbergExchange1D(int site, const std::vector<double> &mxTerms,
                                                               const std::vector<double> &myTerms, const std::vector<double> &mzTerms);

    template <typename T>
    CommonStructures::Vector3D _calculateDMI(int dimension, int site, const std::vector<T> &mxTerms,
                                                       const std::vector<T> &myTerms, const std::vector<T> &mzTerms);

    CommonStructures::Vector3D _calculateDMI1D(int site, const std::vector<double> &mxTerms,
                                                     const std::vector<double> &myTerms, const std::vector<double> &mzTerms);

    template <typename T>
    void addExchangeField(int site, const CommonStructures::Vector3D &exchangeField, std::vector<T> &exchangeXOut,
                                     std::vector<T> &exchangeYOut, std::vector<T> &exchangeZOut);

    bool _hasOscillatingZeeman(int site);

public:
    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<double> &exchangeXOut,
                                std::vector<double> &exchangeYOut, std::vector<double> &exchangeZOut,
                                const CommonStructures::Parallelisations &parallelisationMode);

    void calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms, std::vector<std::atomic<double>> &exchangeXOut,
                                std::vector<std::atomic<double>> &exchangeYOut,
                                std::vector<std::atomic<double>> &exchangeZOut, const CommonStructures::Parallelisations &parallelisationMode);

    CommonStructures::Vector3D calculateExchangeField( int dimension, int site,
                                                  const std::vector<double> &mxTerms,
                                                  const std::vector<double> &myTerms,
                                                  const std::vector<double> &mzTerms);
};


#endif //SPINCHAINS_EXCHANGEFIELD_H
