//
// Created by Cameron McEleney on 01/12/2023.
//

// Corresponding header
#include "../include/ExchangeField.h"

ExchangeField::ExchangeField( SimulationParameters *sharedSimParams,
                              SimulationStates *sharedSimStates,
                              SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {

}

// 1D public
void ExchangeField::calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                           const std::vector<double> &mzTerms, std::vector<double> &exchangeXOut,
                                           std::vector<double> &exchangeYOut, std::vector<double> &exchangeZOut,
                                           const CommonStructures::Parallelisations &parallelisationMode) {
    calculateExchangeField(1, mxTerms, myTerms, mzTerms, exchangeXOut, exchangeYOut, exchangeZOut, parallelisationMode);
}

void ExchangeField::calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                           const std::vector<double> &mzTerms, std::vector<std::atomic<double>> &exchangeXOut,
                                           std::vector<std::atomic<double>> &exchangeYOut,
                                           std::vector<std::atomic<double>> &exchangeZOut,
                                           const CommonStructures::Parallelisations &parallelisationMode) {
    calculateExchangeField(1, mxTerms, myTerms, mzTerms, exchangeXOut, exchangeYOut, exchangeZOut, parallelisationMode);
}

template <typename T>
void ExchangeField::calculateExchangeField( int dimension, const std::vector<double> &mxTerms,
                                            const std::vector<double> &myTerms, const std::vector<double> &mzTerms,
                                            std::vector<T> &exchangeXOut, std::vector<T> &exchangeYOut,
                                            std::vector<T> &exchangeZOut,
                                            CommonStructures::Parallelisations parallelFlag ) {

    // Note: If additional multithreading options are added, the method signature will need to be updated

    if (parallelFlag == CommonStructures::Parallelisations::Multithreaded) {
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins + 1),
            [&](const tbb::blocked_range<int> &range) {
                for (int site = range.begin(); site < range.end(); site++) {
                    auto tempExchangeLocal = calculateExchangeField(dimension, site, mxTerms, myTerms, mzTerms);
                    addExchangeField(site, tempExchangeLocal, exchangeXOut, exchangeYOut, exchangeZOut);
                }
            }, tbb::auto_partitioner());
    }
    else if (parallelFlag == CommonStructures::Parallelisations::Sequential) {
        #pragma unroll
        for (int site = 1; site <= _simParams->systemTotalSpins; site++) {
            auto tempExchangeLocal = calculateExchangeField(dimension, site, mxTerms, myTerms, mzTerms);
            addExchangeField(site, tempExchangeLocal, exchangeXOut, exchangeYOut, exchangeZOut);
        }
    }
    else { throw std::invalid_argument("ExchangeField::calculateExchangeField: Invalid parallelisation flag"); }
}
// One function for the user to access for the 1D case
CommonStructures::Vector3D ExchangeField::calculateExchangeField( int dimension, int site,
                                                             const std::vector<double> &mxTerms,
                                                             const std::vector<double> &myTerms,
                                                             const std::vector<double> &mzTerms) {
    // High-level function that can check the various Boolean flags of the (direct) exchange field components
    CommonStructures::Vector3D heisenbergExchangeTerms{0.0, 0.0, 0.0};

    if (dimension == 1) {
        auto tempExchangeLocal = _calculateHeisenbergExchange1D(site, mxTerms, myTerms, mzTerms);
        heisenbergExchangeTerms.x() += tempExchangeLocal.x();
        heisenbergExchangeTerms.y() += tempExchangeLocal.y();
        heisenbergExchangeTerms.z() += tempExchangeLocal.z();
    }
    else {
        std::cout << "ExchangeField::calculateExchangeField: Dimensionality beyond 1 not supported" << std::endl;
    }

    if (_simFlags->hasDMI) {
        auto tempDMILocal = _calculateDMI(dimension, site, mxTerms, myTerms, mzTerms);
        heisenbergExchangeTerms.x() += tempDMILocal.x();
        heisenbergExchangeTerms.y() += tempDMILocal.y();
        heisenbergExchangeTerms.z() += tempDMILocal.z();
    }

    return heisenbergExchangeTerms;
}

CommonStructures::Vector3D ExchangeField::_calculateHeisenbergExchange1D(int site, const std::vector<double> &mxTerms,
                                                               const std::vector<double> &myTerms, const std::vector<double> &mzTerms) {

    // Better to pre-compute the indexes and recycle the values each time
    CommonStructures::VectorND<int, 2> siteLocation{site - 1, site + 1};

    // Use references to avoid creating additional memory; want separate variables to improve readability
    auto &exchangeLhs = _simStates->exchangeVec[siteLocation.x()];
    auto &exchangeRhs = _simStates->exchangeVec[site];

    CommonStructures::Vector3D heisenbergExchangeTerms{0.0, 0.0, 0.0};

    heisenbergExchangeTerms.x() = exchangeLhs * mxTerms[siteLocation.x()] + exchangeRhs * mxTerms[siteLocation.y()];
    heisenbergExchangeTerms.y() = exchangeLhs * myTerms[siteLocation.x()] + exchangeRhs * myTerms[siteLocation.y()];
    heisenbergExchangeTerms.z() = exchangeLhs * mzTerms[siteLocation.x()] + exchangeRhs * mzTerms[siteLocation.y()];

    if (!_simFlags->isFerromagnetic) {
        heisenbergExchangeTerms.x() *= -1.0;
        heisenbergExchangeTerms.y() *= -1.0;
        heisenbergExchangeTerms.z() *= -1.0;
        heisenbergExchangeTerms.z() += (mzTerms[site] > 0) ? _simParams->anisotropyField : - _simParams->anisotropyField;
    }

    return heisenbergExchangeTerms;
}

template <typename T>
CommonStructures::Vector3D ExchangeField::_calculateDMI(int dimension, int site, const std::vector<T> &mxTerms,
                                                   const std::vector<T> &myTerms, const std::vector<T> &mzTerms) {

    CommonStructures::Vector3D dmiExchangeTerms{0.0, 0.0, 0.0};  // Container for output
    if (dimension == 1) { dmiExchangeTerms = _calculateDMI1D(site, mxTerms, myTerms, mzTerms); }
    else { throw std::invalid_argument("ExchangeField::_calculateDMI: Dimensionality beyond 1 not supported"); }

    // Track scaling due to site's location within the gradient region
    double scalingFactor = 1.0;  // Currently `_simStates->dmiGradientMap` is a `double` type

    // Use map to check if given site is within stated gradient region
    auto iter = _simStates->dmiGradientMap.find(site);
    if ( iter != _simStates->dmiGradientMap.end()) { scalingFactor = iter->second; }
    else {
        // Site isn't in gradient region, so it's here we must check further criteria
        if ( _simFlags->shouldRestrictDmiToWithinGradientRegion ) { scalingFactor = 0.0; }
    }

    // If scaling factor is unchanged for current site after gradient map checks, then we have an early return opportunity
    if (scalingFactor == 1.0) { return dmiExchangeTerms; }

    #pragma unroll 3
    for (int i = 0; i < 3; i++) { dmiExchangeTerms[i] *= scalingFactor; }

    return dmiExchangeTerms;
}

CommonStructures::Vector3D ExchangeField::_calculateDMI1D(int site, const std::vector<double> &mxTerms,
                                                     const std::vector<double> &myTerms, const std::vector<double> &mzTerms) {
    
    CommonStructures::Vector3D siteLhs = {mxTerms[site - 1], myTerms[site - 1], mzTerms[site - 1]};
    CommonStructures::Vector3D siteRhs = {mxTerms[site + 1], myTerms[site + 1], mzTerms[site + 1]};

    CommonStructures::Vector3D dmiEquation = {
         _simParams->dmiVector.y() * (siteRhs.z() - siteLhs.z()) - _simParams->dmiVector.z() * (siteRhs.y() - siteLhs.y()),
        -_simParams->dmiVector.x() * (siteRhs.z() - siteLhs.z()) + _simParams->dmiVector.z() * (siteRhs.x() - siteLhs.x()),
         _simParams->dmiVector.x() * (siteRhs.y() - siteLhs.y()) - _simParams->dmiVector.y() * (siteRhs.x() - siteLhs.x())
    };


    /* Equivalent to
     * return {_simParams->dmiConstant * (myTerms[currentSite + 1] - myTerms[currentSite - 1]),
     *         -1.0 * _simParams->dmiConstant * (mxTerms[currentSite + 1] - mxTerms[currentSite - 1]),
     *         0.0}
     */

    return dmiEquation;
}

template <typename T>
void ExchangeField::addExchangeField(int site, const CommonStructures::Vector3D &exchangeField, std::vector<T> &exchangeXOut,
                                     std::vector<T> &exchangeYOut, std::vector<T> &exchangeZOut) {
    exchangeXOut[site] += exchangeField.x();
    exchangeYOut[site] += exchangeField.y();
    exchangeZOut[site] += exchangeField.z();
}

// Specialization for std::atomic<double>
template <>
void ExchangeField::addExchangeField<std::atomic<double>>(int site, const CommonStructures::Vector3D &exchangeField,
                                                          std::vector<std::atomic<double>> &exchangeXOut,
                                                          std::vector<std::atomic<double>> &exchangeYOut,
                                                          std::vector<std::atomic<double>> &exchangeZOut) {
    exchangeXOut[site].fetch_add(exchangeField.x(), std::memory_order_relaxed);
    exchangeYOut[site].fetch_add(exchangeField.y(), std::memory_order_relaxed);
    exchangeZOut[site].fetch_add(exchangeField.z(), std::memory_order_relaxed);
}

// Cross product method
CommonStructures::Vector3D ExchangeField::_crossProduct(const CommonStructures::Vector3D &iSite,
                                                   const CommonStructures::Vector3D &jSite) {

    return { iSite.y() * jSite.z() - iSite.z() * jSite.y(),
            -iSite.x() * jSite.z() + iSite.z() * jSite.x(),
             iSite.x() * jSite.y() - iSite.y() * jSite.x()
    };
}

bool ExchangeField::_hasOscillatingZeeman(int site) {
    if (_simFlags->shouldDriveDiscreteSites) {
        for (const int &discreteSite: _simStates->discreteDrivenSites) {
            if (site == discreteSite) { return true; }
        }
    }

    if (site >= _simParams->drivingRegion.lhsSite && site <= _simParams->drivingRegion.rhsSite) { return true; }

    // If no condition is met, then the site is not driven
    return false;
}

/*
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
            [&]( const tbb::blocked_range<int> &tbbRange ) {
                for ( int site = tbbRange.begin(); site < tbbRange.end(); ++site ) { // for some reason site++ doesn't work
                    CommonStructures::Vector3D tempResultsExchangeLocal = _calculateExchangeField1D(site, mxTerms, myTerms, mzTerms, shouldUseTBB);

                    // Use of 'auto' allows for tertiary operator to be used; equivalent to declaring array and then initialising within an IF/ELSE structure
                    auto tempResultsDMILocal = _simFlags->hasDMI ? _calculateDMI1D(site, mxTerms, myTerms, mzTerms, shouldUseTBB)
                                                                 : CommonStructures::Vector3D{0.0, 0.0, 0.0};

                    tbb::mutex::scoped_lock lock;  // Needed to prevent race condition
                    exchangeXOut[site] += tempResultsExchangeLocal.x() + tempResultsDMILocal.x();
                    exchangeYOut[site] += tempResultsExchangeLocal.y() + tempResultsDMILocal.y();
                    exchangeZOut[site] += tempResultsExchangeLocal.z() + tempResultsDMILocal.z();
                }
            }, tbb::auto_partitioner());


        // Thread-local storage for each component
        tbb::combinable<std::vector<double>> localFieldX([&]{ return std::vector<double>(_simParams->systemTotalSpins + 2, 0.0); });
        tbb::combinable<std::vector<double>> localFieldY([&]{ return std::vector<double>(_simParams->systemTotalSpins + 2, 0.0); });
        tbb::combinable<std::vector<double>> localFieldZ([&]{ return std::vector<double>(_simParams->systemTotalSpins + 2, 0.0); });

        // Parallel computation
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
            [&](const tbb::blocked_range<int>& range) {
                auto& localX = localFieldX.local();
                auto& localY = localFieldY.local();
                auto& localZ = localFieldZ.local();

                for (int i = range.begin(); i <= range.end(); i++) {
                    CommonStructures::Vector3D fieldContribution = _calculateExchangeField1D(i, mxTerms, myTerms, mzTerms, shouldUseTBB);
                    localX[i] += fieldContribution.x();
                    localY[i] += fieldContribution.y();
                    localZ[i] += fieldContribution.z();
                }
            }, tbb::auto_partitioner());

        // Combining the results directly into the output fields
        localFieldX.combine_each([&](const std::vector<double>& v) {
            for (size_t i = 1; i <= _simParams->systemTotalSpins; i++) exchangeXOut[i] += v[i];
        });
        localFieldY.combine_each([&](const std::vector<double>& v) {
            for (size_t i = 1; i <= _simParams->systemTotalSpins; i++) exchangeYOut[i] += v[i];
        });
        localFieldZ.combine_each([&](const std::vector<double>& v) {
            for (size_t i = 1; i <= _simParams->systemTotalSpins; i++) exchangeZOut[i] += v[i];
        });
         */