//
// Created by Cameron McEleney on 01/12/2023.
//

#ifndef SPINCHAINS_BIASFIELDS_H
#define SPINCHAINS_BIASFIELDS_H


// C++ Standard Library

// C++ Third Party Libraries
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

// C++ User Libraries (General)
#include "CommonLibs.h"
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"

/*
 * This class should calculate all of the bias fields required by the system. As there might be many of these,
 * but they are often simple, it is worth while grouping them together as a 'single' component of the system.
 *
 * Bias fields are external fields that influence the system behaviour. They are also often encountered under the name of
 * Zeeman fields.
 *
 * Examples of Zeeman fields, which are all bias field, and that could be included are:
 *      - Static Field. Also called 'external static'. Often labelled 'H_0', and used to establish a preferred magnetisaton direction.
 *              - Can be uniform or non-uniform, like \mu_0 * H_0 = 0.1 [T] in macdo2022breaking.
 *
 *      - Oscillating Field. Also called 'driving field' and 'external bias'.
 *              - Often time-dependent, like h(x, t) = 3[mT] * sin(2\pi * f_d * t) [T] in macedo2022breaking.
 *              - Can be uniform or non-uniform
 *
 *      - Gradient Field
 *              - Used to simulate effects of magnetic field gradients, and is spatially varying. Used in magnetic recording,
 *              such as HAMR, and in MRIs.
 *
 *      - Pulsed Field.
 *              - Short, often ultrafast, bursts of a magnetic field used to manipulate the magnetisation like
 *                in magnetic memory (MRAM) applications or to manipulate transition between energy levels.
 *
 * Examples of bias fields, which are not Zeeman fields, and that could be included are:
 *      - Random Fields
 *              - Simulate (external) thermal fluctuations or noise that enter the system from the environment. Different
 *                to effects from internal thermal effects!
 */

class BiasFields {
private:
    SimulationParameters *_simParams; // Non-owning pointer to SimulationParameters
    SimulationStates *_simStates;
    SimulationFlags *_simFlags;

private:
    std::array<double, 3>
    _calculateBiasField1D( const int &currentSite, const int &currentLayer, const double &currentTime );

    std::array<double, 3>
    _calculateBiasField1D( const int &currentSite, const int &currentLayer, const double &currentTime,
                           const double &mzTermAtSite, const bool &shouldUseTBB );

    inline bool _hasOscillatingZeeman( const int &site );

public:
    explicit BiasFields( SimulationParameters *sharedSimParams,
                         SimulationStates *sharedSimStates,
                         SimulationFlags *sharedSimFlags );

    ~BiasFields() = default;

public:
    void calculateOneDimension( const int &currentLayer, const double &currentTime, const std::vector<double> &mzTermsIn,
                                std::vector<double> &biasFieldXOut, std::vector<double> &biasFieldYOut,
                                std::vector<double> &biasFieldZOut );

    void calculateOneDimension( const int &currentLayer, const double &currentTime,
                                const std::vector<double> &mzTermsIn, std::vector<std::atomic<double>> &biasFieldXOut,
                                std::vector<std::atomic<double>> &biasFieldYOut, std::vector<std::atomic<double>> &biasFieldZOut,
                                const bool &shouldUseTBB );

    void calculateOneDimension( const int &currentLayer, const double &currentTime,
                                const std::vector<double> &mzTermsIn, std::vector<double> &biasFieldXOut,
                                std::vector<double> &biasFieldYOut, std::vector<double> &biasFieldZOut,
                                const bool &shouldUseTBB );
};


#endif //SPINCHAINS_BIASFIELDS_H
