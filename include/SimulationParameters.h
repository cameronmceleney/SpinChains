//
// Created by Cameron Aidan McEleney on 19/11/2023.
//

#ifndef SPINCHAINS_SIMULATIONPARAMETERS_H
#define SPINCHAINS_SIMULATIONPARAMETERS_H

// C++ Standard Library
#include "../libs/CommonStructures.h"

struct HeisenbergExchange {
    // Member variables
    float uniformStrength = 0.0;
    std::optional<float> minimumStrength;
    std::optional<float> maximumStrength;

    CommonStructures::Unit units{CommonStructures::Unit::Tesla};
    bool isExchangeUniform = true;

    // Default constructor
    HeisenbergExchange() = default;

    // Constructor for uniform exchange
    HeisenbergExchange(float value, CommonStructures::Unit unit)
        : uniformStrength(value), units(unit), isExchangeUniform(true) {
            minimumStrength.reset();
            maximumStrength.reset();
        }

    // Constructor for non-uniform exchange
    HeisenbergExchange(float minStrength, float maxStrength, CommonStructures::Unit unit)
        : minimumStrength(minStrength), maximumStrength(maxStrength), units(unit), isExchangeUniform(false) {}

    // Method to set uniform exchange (replaces GlobalVariables.cpp)
    void setUniformExchange(float value, CommonStructures::Unit unit) {
        uniformStrength = value;
        units = unit;
        isExchangeUniform = true;
        minimumStrength.reset();
        maximumStrength.reset();

        std::cout << "Uniform exchange interaction set." << std::endl;
    }

    // Method to set non-uniform exchange (replaces GlobalVariables.cpp)
    void setNonUniformExchange(float minStrength, float maxStrength, CommonStructures::Unit unit) {
        minimumStrength = minStrength;
        maximumStrength = maxStrength;
        units = unit;
        isExchangeUniform = false;
        uniformStrength = 0.0;

        std::cout << "Non-uniform exchange interaction set." << std::endl;
    }

    // Getters to replace contents in GlobalVariables.cpp
    [[nodiscard]] std::pair<float, CommonStructures::Unit> getUniformExchange() const {
        if (!isExchangeUniform) {
            throw std::logic_error("Exchange interaction is not uniform");
        }
        return {uniformStrength, units};
    }

    [[nodiscard]] std::tuple<float, float, CommonStructures::Unit> getNonUniformExchange() const {
        if (isExchangeUniform) {
            throw std::logic_error("Exchange interaction is uniform");
        }
        if (!minimumStrength.has_value() || !maximumStrength.has_value()) {
            throw std::logic_error("Non-uniform exchange interaction values are not set");
        }
        return {*minimumStrength, *maximumStrength, units};
    }
};

struct GilbertDamping {
    double factor;  // Gilbert damping factor for main chain.
    double innerABC;  // The lower Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the chain meets the ABC.
    double outerABC;  // The upper Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the ABC meets the pinned sites.

    GilbertDamping() = default;
};

struct MappedRegion {
    int lhsSite = -1;
    int rhsSite = -1;
    int peakWidth;
    int gradientWidth;
    int regionWidth = -1;  // This should be automatically calculated
    int offsetWidth;

    double peakValue;
    double gradientValue;
};

class SimulationParameters {
    // todo separate this out into several compositions (flags, data structures, etc)
public:
    constexpr static double    PERM_FREESPACE = 1.25663706212e-6;         // Permeability of free space [H/m] (mu_0)
    constexpr static double    BOLTZMANN_CONSTANT = 1.380649e-23;         // Boltzmann Constant [m^{2} kg s^{-2} K^{-1}].
    double PERMEABILITY_IRON = 2.22; // cobalt = 1.72, iron = 2.22, nickel = 0.6;

    // Sites to be printed if shouldPrintDiscreteSites is TRUE.


    double              ambientTemperature;
    double              anisotropyField;
    CommonStructures::Vector3D  staticZeemanStrength{0.0, 0.0, 0.0};

    double              drivingAngFreq;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
    double              drivingFreq;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)

    double              oscillatingZeemanStrength;                // Driving field amplitude [T] (caution: papers often give in [mT]).
    int                 forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
    double              dipoleConstant;                           // Scaling factor which is constant across dipolar interaction calculations.
    CommonStructures::Vector3D  dmiVector{0.0, 0.0, 0.0};
    double              dmiConstant;
    double              exchangeStiffness;

    double              gyroMagConst;                             // Gyromagnetic ratio of an electron [GHz/T].
    double              latticeConstant;

    double              spinPolarisation;                         // Spin polarisation of the spin current.
    double              spinTransferEfficiency;                   // Spin transfer efficiency of the spin current.

    int                 iterationEnd;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    int                 iterationStart = 0;                        // The iteration step that the program will begin at (Default: 0.0)
    double              risingTimeStartAtIteration;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the maxSimTime.

    double              risingTimeEndAtIteration;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the maxSimTime.
    double              largestMNorm = 1e-50;                     // Computes sqrt(_mx**2 + _my**2 + _mz**2) for each site at each moment to track any abnormalities. Initialises to be arbitrarily small
    double              maxSimTime;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [mxInit + myInit + mzInit]  CANNOT sum to greater than 1.0
    int                 numNeighbours;
    int                 numberOfDataPoints;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    int                 numberOfSpinPairs;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    int                 numSpinsInABC;                           // Number of spins in the damped regions (previously called _numGilbert).

    int                 numSpinsInChain;                          // The number of spin sites in the spin chain to be simulated.
    int                 systemTotalSpins;                         // The total number of spins in the system (chain plus ABCs).
    double              satMag;                                   // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)

    double              risingTimePeriod;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    double              risingTimeInitialMagnitude;                 // Initial strength of the shockwave before risingTimeScalingFactor occurs. (Default: = oscillatingZeemanStrength)
    double              risingTimeMaximum;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    double              risingTimeScalingFactor;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving

    double              shockwaveStepsize;                        // Size of incremental increase in shockwave amplitude.
    double              stepsize;                                 // Stepsize between values
    double              stepsizeHalf;                             // Separately defined to avoid repeated unnecessary calculations inside loops
    int                 numLayers;
    double              recordingInterval;                        // Time between each data point recorded.
    int                 layerOfInterest;

    double exchangeEnergyMin;
    double exchangeEnergyMax;

    double totalTime = 0;

    //TODO. Turn the below lines into structs (much easier to access data)
    int numSpinsDRPeak, numSpinsDRGradient;

public:
    HeisenbergExchange exchangeInteraction;
    GilbertDamping gilbertDamping;
    MappedRegion dmiRegion;
    MappedRegion drivingRegion;
    MappedRegion dampingRegion;
};

#endif //SPINCHAINS_SIMULATIONPARAMETERS_H
