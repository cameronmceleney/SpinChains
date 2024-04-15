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

        exchangeXOut[i] = tempResultsContainerHeisenberg[0] + tempResultsContainerDMI[0];
        exchangeYOut[i] = tempResultsContainerHeisenberg[1] + tempResultsContainerDMI[1];
        exchangeZOut[i] = tempResultsContainerHeisenberg[2] + tempResultsContainerDMI[2];
    }
}


void ExchangeField::calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                           const std::vector<double> &mzTerms, std::vector<double> &effectiveFieldX,
                                           std::vector<double> &effectiveFieldY, std::vector<double> &effectiveFieldZ,
                                           const int &selectThisFunction) {
    /*
     * mxTerms/myTerms/mzTerms take us all the back to the class-wide containers for the m-components
     * from SolversImplementation.cpp. As the actual containers will change depending on RK2 stage, we can't access
     * the class container via this utility method directly.
     *
     * effectiveFieldX/effectiveFieldY/effectiveFieldZ come from class-wide containers from
     * SolversImplementation.cpp. Can't allow direct access though in case of higher level divide-and-conquer;
     * calculateOneDimension is just a component so it should only work on what it is passed.
     *
     * We pass-by-reference to follow general C++ style.
     */

    // Value of each element should be 3 as only 3 reads are required; we load mx/my/mz as a set so we only need to count
    // the number of sets. This turns 3 vectors of counters into a single vector counter
    //std::vector<std::atomic<int16_t>> syncVec(_simParams->systemTotalSpins + 2);
    //for (auto& atomicInt : syncVec) {
    //    atomicInt.store(3, std::memory_order_relaxed);
    //}

    // Adjust the boundary conditions: sites at the edges have only 2 reads
    //syncVec[1].store(2, std::memory_order_relaxed);
    //syncVec[_simParams->systemTotalSpins + 1].store(2, std::memory_order_relaxed);

    //std::atomic<int> active_threads(0);

    //tbb::global_control c( tbb::global_control::max_allowed_parallelism, 1 );
    tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins + 1),
        [&](const tbb::blocked_range<int> &range) {
            for (int site = range.begin(); site < range.end(); site++) {
                //active_threads.fetch_add(1, std::memory_order_relaxed);
                // Single temporary container for this thread's results
                std::array<double, 3> tempExchangeLocal{0.0, 0.0, 0.0};

                // Returns the exchangeField updates in tempExchangeLocal
                _calculateExchangeField1D(site, tempExchangeLocal, mxTerms, myTerms, mzTerms);

                if (_simFlags->hasDMI)
                {
                    // Returns the dmiField updates in tempExchangeLocal
                    _calculateDMI1D(site, tempExchangeLocal,
                                    mxTerms, myTerms, mzTerms);
                }

                /*
                if (_simFlags->hasDemag1DThinFilm) {
                    // Reduces total operations by only summing when DMI is present
                    _calculateDemagSimple(site, tempExchangeLocal, mxTerms, myTerms, mzTerms);
                }
                 */
                // Decrement the counters for the adjacent sites, signaling this thread's completion of their usage.
                //syncVec[site-1].fetch_sub(1, std::memory_order_release);
                //syncVec[site+1].fetch_sub(1, std::memory_order_release);

                // Decrement the counter for the current site and check if this thread is the last to need its data.
                //int countsRemaining = syncVec[site].fetch_sub(1, std::memory_order_acq_rel) - 1;

                // If this thread is not the last one (remaining > 0), it needs to wait.
                //if (countsRemaining > 0) {
                //    auto start = std::chrono::steady_clock::now();
                //    auto now = std::chrono::steady_clock::now();
                //    while (syncVec[site].load(std::memory_order_acquire) != 0) {
                //         now = std::chrono::steady_clock::now();
                //        if (std::chrono::duration_cast<std::chrono::seconds>(now - start).count() > 5) {
                //            // Timeout exceeded, print the values within syncVec
                //            std::cout << "Timeout exceeded, syncVec values are: \n";
                //            for (size_t i = 0; i < syncVec.size(); ++i) {
                //                std::cout << "syncVec[" << i << "]: " << syncVec[i].load(std::memory_order_acquire) << "\n";
                //            }
                //            std::exit(active_threads); // Break out of the loop or handle the situation as needed
                //        }
                //    }
                //}

                // Safe to update exchangeVec. Remember that this is the class-wide container!
                effectiveFieldX[site] = tempExchangeLocal[0];
                effectiveFieldY[site] = tempExchangeLocal[1];
                effectiveFieldZ[site] = tempExchangeLocal[2];
                // Reset syncVec[site] here if its defined out with the function scope
                //active_threads.fetch_sub(1, std::memory_order_relaxed);

            }
    }, tbb::auto_partitioner());
    //bufferEffectiveFieldX = std::move(effectiveFieldX);
    //bufferEffectiveFieldY = std::move(effectiveFieldY);
    //bufferEffectiveFieldZ = std::move(effectiveFieldZ);

}

void ExchangeField::calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                           const std::vector<double> &mzTerms, std::vector<double> &exchangeXOut,
                                           std::vector<double> &exchangeYOut, std::vector<double> &exchangeZOut,
                                           const bool &shouldUseTBB ) {
    // Doesn't work right now


    if ( shouldUseTBB ) {
        /*
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
            [&]( const tbb::blocked_range<int> &tbbRange ) {
                for ( int site = tbbRange.begin(); site < tbbRange.end(); ++site ) { // for some reason site++ doesn't work
                    std::array<double, 3> tempResultsExchangeLocal = _calculateExchangeField1D(site, mxTerms, myTerms, mzTerms, shouldUseTBB);

                    // Use of 'auto' allows for tertiary operator to be used; equivalent to declaring array and then initialising within an IF/ELSE structure
                    auto tempResultsDMILocal = _simFlags->hasDMI ? _calculateDMI1D(site, mxTerms, myTerms, mzTerms, shouldUseTBB)
                                                                 : std::array<double, 3>{0.0, 0.0, 0.0};

                    tbb::mutex::scoped_lock lock;  // Needed to prevent race condition
                    exchangeXOut[site] += tempResultsExchangeLocal[0] + tempResultsDMILocal[0];
                    exchangeYOut[site] += tempResultsExchangeLocal[1] + tempResultsDMILocal[1];
                    exchangeZOut[site] += tempResultsExchangeLocal[2] + tempResultsDMILocal[2];
                }
            }, tbb::auto_partitioner());
        */

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
                    std::array<double, 3> fieldContribution = _calculateExchangeField1D(i, mxTerms, myTerms, mzTerms, shouldUseTBB);
                    localX[i] += fieldContribution[0];
                    localY[i] += fieldContribution[1];
                    localZ[i] += fieldContribution[2];
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
    } else {
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

    /* for some reason this commented code doesn't work
     *     if (shouldUseTBB) {
        if ( _simFlags->isFerromagnetic ) {
            return {
                _simStates->exchangeVec[currentSite - 1] * mxTerms[currentSite - 1] + _simStates->exchangeVec[currentSite] * mxTerms[currentSite + 1],
                _simStates->exchangeVec[currentSite - 1] * myTerms[currentSite - 1] + _simStates->exchangeVec[currentSite] * myTerms[currentSite + 1],
                _simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1] + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1]
            };
        } else if ( !_simFlags->isFerromagnetic ) {

            double mzAnisotropyTerm; // Written this way to be explicit for the anisotropy term
            if ( mzTerms[currentSite] > 0 ) { mzAnisotropyTerm = _simParams->anisotropyField; }
            else if ( mzTerms[currentSite] < 0 ) {mzAnisotropyTerm =  -1.0 * _simParams->anisotropyField; }

            return {
                -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTerms[currentSite - 1] + _simStates->exchangeVec[currentSite] * mxTerms[currentSite + 1]),
                -1.0 * (_simStates->exchangeVec[currentSite - 1] * myTerms[currentSite - 1] + _simStates->exchangeVec[currentSite] * myTerms[currentSite + 1]),
                mzAnisotropyTerm - (_simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1] + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1])
            };
        } else {
            throw std::invalid_argument("_calculateExchangeField1D hasn't got CUDA implementation yet");
        }
    }
     */

    if ( shouldUseTBB ) {
        // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
        double directExchangeXLocal, directExchangeYLocal, directExchangeZLocal;

        // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

        if ( _simFlags->isFerromagnetic ) {
            // hx terms
            directExchangeXLocal = _simStates->exchangeVec[currentSite - 1] * mxTerms[currentSite - 1]
                        + _simStates->exchangeVec[currentSite] * mxTerms[currentSite + 1];

            // hy terms
            directExchangeYLocal = _simStates->exchangeVec[currentSite - 1] * myTerms[currentSite - 1]
                          + _simStates->exchangeVec[currentSite] * myTerms[currentSite + 1];

            // hz terms
            directExchangeZLocal = _simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1]
                          + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1];

        } else {
            // hx terms
            directExchangeXLocal = -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTerms[currentSite - 1]
                        + _simStates->exchangeVec[currentSite] * mxTerms[currentSite + 1]);


            // hy terms
            directExchangeYLocal = -1.0 * (_simStates->exchangeVec[currentSite - 1] * myTerms[currentSite - 1]
                                  + _simStates->exchangeVec[currentSite] * myTerms[currentSite + 1]);

            // hz terms
            if ( mzTerms[currentSite] > 0 )
                directExchangeZLocal = _simParams->anisotropyField -
                              (_simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1]
                               + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1]);
            else if ( mzTerms[currentSite] < 0 )
                directExchangeZLocal =  -1.0 * _simParams->anisotropyField -
                              (_simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1]
                               + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1]);
        }

        return {directExchangeXLocal, directExchangeYLocal, directExchangeZLocal};
    } else
        throw std::invalid_argument("_calculateExchangeField1D hasn't got CUDA implementation yet");
}

std::array<double, 3>&
ExchangeField::_calculateExchangeField1D( const int &currentSite, std::array<double, 3> &directExchangeLocal,
                                          const std::vector<double> &mxTerms,
                                          const std::vector<double> &myTerms,
                                          const std::vector<double> &mzTerms ) {

    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double directExchangeXLocal, directExchangeYLocal, directExchangeZLocal;

    // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

    if ( _simFlags->isFerromagnetic ) {
        // hx terms
        directExchangeLocal[0] = _simStates->exchangeVec[currentSite - 1] * mxTerms[currentSite - 1]
                    + _simStates->exchangeVec[currentSite] * mxTerms[currentSite + 1];

        // hy terms
        directExchangeLocal[1] = _simStates->exchangeVec[currentSite - 1] * myTerms[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * myTerms[currentSite + 1];

        // hz terms
        directExchangeLocal[2] = _simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1]
                      + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1];

    } else {
        // hx terms
        directExchangeLocal[0] = -1.0 * (_simStates->exchangeVec[currentSite - 1] * mxTerms[currentSite - 1]
                    + _simStates->exchangeVec[currentSite] * mxTerms[currentSite + 1]);


        // hy terms
        directExchangeLocal[1] = -1.0 * (_simStates->exchangeVec[currentSite - 1] * myTerms[currentSite - 1]
                              + _simStates->exchangeVec[currentSite] * myTerms[currentSite + 1]);

        // hz terms
        if ( mzTerms[currentSite] > 0 )
            directExchangeLocal[2] = _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1]);
        else if ( mzTerms[currentSite] < 0 )
            directExchangeLocal[2] =  -1.0 * _simParams->anisotropyField -
                          (_simStates->exchangeVec[currentSite - 1] * mzTerms[currentSite - 1]
                           + _simStates->exchangeVec[currentSite] * mzTerms[currentSite + 1]);
    }

    return directExchangeLocal;
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
        return {_simParams->dmiConstant * (myTerms[currentSite + 1] - myTerms[currentSite - 1]),
                 -1.0 * _simParams->dmiConstant * (mxTerms[currentSite + 1] - mxTerms[currentSite - 1]),
                 0.0};
    } else {
        throw std::invalid_argument("_calculateDMIField1D hasn't got CUDA implementation yet");
    }
}

std::array<double, 3>&
ExchangeField::_calculateDMI1D( const int &currentSite, std::array<double, 3>& localDmiTerms,
                                const std::vector<double> &mxTerms,
                                const std::vector<double> &myTerms,
                                const std::vector<double> &mzTerms ) {

    /* TODO. This is a temp polymorphic version of _calculateDMIField1D that is threadsafe. Needs refinement
     * Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
     *
     * Quickly finds solution to H_DMI = D_{i-1, i} \cdot (m_{i-1} \times m_{i}) + D_{i, i+1} \cdot (m_{i} \times m_{i+1}
     * by exploiting how only D_z is non-zero in this 1D system
     *
     * Uses Eq. 3 from https://doi.org/10.1103/PhysRevB.107.224421 to explicitly write the return statements for the case
     * where DMI only involves two nearest neighbours (NN).
     */
    // Initialisation serves as ELSE: No gradient for DMI, or, site not in map and not restricting dmi to
    // gradient region; linear DMI scaling throughout system
     double scaledDmiConstant = _simParams->dmiConstant;

    if ( _simFlags->hasGradientRegionForDmi )
    {
        auto it = _simStates->dmiGradientMap.find(currentSite);

        if ( it != _simStates->dmiGradientMap.end())
        {
            // Found the site in the map;
            scaledDmiConstant *= it->second;
        }
        else if ( _simFlags->shouldRestrictDmiToWithinGradientRegion )
        {
            // Site not in the map
            scaledDmiConstant *= 0.0;
        }
    }

    localDmiTerms[0] += scaledDmiConstant  * (myTerms[currentSite + 1] - myTerms[currentSite - 1]);
    localDmiTerms[1] += -1.0 * scaledDmiConstant  * (mxTerms[currentSite + 1] - mxTerms[currentSite - 1]);
    // localDmiTerms[2] += 0.0;

    return localDmiTerms;
}

std::array<double, 3>
ExchangeField::_calculateDemagSimple( const int &currentSite, const std::vector<double> &mxTerms,
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
        return {-1.005096 * myTerms[currentSite] * mzTerms[currentSite],
                0,
                -1.005096 * myTerms[currentSite] * mxTerms[currentSite]};
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

bool ExchangeField::_hasOscillatingZeeman( const int &site ) {
    if ( _simFlags->shouldDriveDiscreteSites ) {
        for ( const int &discreteSite: _simStates->discreteDrivenSites )
            if ( site == discreteSite ) { return true; }
    }

    if ( site >= _simParams->drivingRegionLhs && site <= _simParams->drivingRegionRhs ) { return true; }

    // If no condition is met, then the site is not driven
    return false;
}
