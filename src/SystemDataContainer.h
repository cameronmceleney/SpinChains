//
// Created by Cameron Aidan McEleney on 19/11/2023.
//

#include <iomanip>
#include <list>

#include "../src/CommonLibs.h"

#ifndef SPINCHAINS_SYSTEMDATACONTAINER_H
#define SPINCHAINS_SYSTEMDATACONTAINER_H

class SystemDataContainer {
    // todo separate this out into several compositions (flags, data structures, etc)
public:
    double              ambientTemperature;
    double              anisotropyField;
    constexpr static double    BOLTZMANN_CONSTANT = 1.380649e-23;         // Boltzmann Constant [m^{2} kg s^{-2} K^{-1}].

    double              drivingAngFreq;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
    double              drivingFreq;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
    int                 drivingRegionLhs;                         // The position of the spin which is leftmost in the driving region.
    int                 drivingRegionRhs;                         // The position of the spin which is rightmost in the driving region.

    int                 drivingRegionWidth;                       // Driving region width.
    double              dynamicBiasField;                         // Driving field amplitude [T] (caution: papers often give in [mT]).
    int                 forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
    double              dipoleConstant;                           // Scaling factor which is constant across dipolar interaction calculations.

    double              gilbertDamping;                           // Gilbert damping factor for main chain.
    double              gilbertABCInner;                          // The lower Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the chain meets the ABC.
    double              gilbertABCOuter;                          // The upper Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the ABC meets the pinned sites.
    double              gyroMagConst;                             // Gyromagnetic ratio of an electron [GHz/T].

    int                 iterationEnd;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    double              iterStartShock;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the maxSimTime.

    double              iterEndShock;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the maxSimTime.
    double              maxSimTime;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [mxInit + myInit + mzInit]  CANNOT sum to greater than 1.0
    int                 numberNeighbours;
    int                 numberOfDataPoints;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    int                 numberOfSpinPairs;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    int                 numSpinsDamped;                           // Number of spins in the damped regions (previously called _numGilbert).

    int                 numSpinsInChain;                          // The number of spin sites in the spin chain to be simulated.
    int                 systemTotalSpins;                         // The total number of spins in the system (chain plus ABCs).
    constexpr static double    PERM_FREESPACE = 1.25663706212e-6;         // Permeability of free space [H/m] (mu_0)
    double              satMag;                                   // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)

    double              shockwaveGradientTime;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    double              shockwaveInitialStrength;                 // Initial strength of the shockwave before shockwaveScaling occurs. (Default: = dynamicBiasField)
    double              shockwaveMax;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    double              shockwaveScaling;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving

    double              shockwaveStepsize;                        // Size of incremental increase in shockwave amplitude.
    double              stepsize;                                 // Stepsize between values
    double              stepsizeHalf;                             // Separately defined to avoid repeated unnecessary calculations inside loops
    std::string         stepsizeString;                           // Object to string conversation for stepsize

    std::string         stopIterString;                           // Object to string conversion for _stopIterVal
    int                 totalLayers;


    // #####################################        Protected Flags        ####################################
    bool                centralDrive;                               // Drive from the centre of the chain if (true)
    bool                driveAllLayers;
    bool                dualDrive;                                  // Drive from both sides of the system
    bool                hasShockwave;                               // Simulation contains a single driving bias field if (false).

    bool                hasStaticDrive;                             // Selects (if true) whether drive has sinusoidal term
    bool                isShockwaveOn;                              // Tests if the conditions to trigger a shockwave have been reached. Not to be altered by the user.
    bool                isShockwaveAtMax;                           // Tests if the shockwave is at its maximum amplitude. Not to be altered by the user.

    bool                lhsDrive;                                   // Drive from the LHS
    bool                rhsDrive;                                   // Drive from the RHS
    bool                printAllData;                               // Saves the m-component(s) of every spin at every iteration. WARNING: leads to huge output files.
    bool                printFixedLines;                            // Saves m-component(s) of every spin at regular intervals. Total save points are set by numberOfDataPoints.
    bool                printFixedSites;                            // Saves a discrete set of m-component(s) at regular intervals governed by numberOfDataPoints.

    bool                shouldDriveCease;                           // Internal flag to indicate if the driving field should cut off at a given time.
    bool                shouldTrackMValues;                         // Monitor the norm of all the m-values; if approx. 1.0 then the error is likely to be massive; discard that dataset.
    bool                useLLG;                                     // Uses the Torque equation components if (false).
    bool                useSLLG;

    bool                useDipolar;
    bool                useZeeman;
    bool                useDemagIntense;                            // todo doesn't work
    bool                useDemagFft;                                // todo doesn't work
    bool                useMultilayer;
    bool                debugFunc;

    bool isFm = GV.GetIsFerromagnetic();
    double exchangeEnergyMin = GV.GetExchangeMinVal();
    double exchangeEnergyMax = GV.GetExchangeMaxVal();

    std::vector<double> gilbertVector = {0};
    std::vector<std::vector<double>> gilbertVectorMulti = {};
    std::vector<double> largestMNormMulti = {1e-50, 1e-50};
    std::vector<double> mx0 = {0};
    std::vector<double> my0 = {0};
    std::vector<double> mz0 = {0};

    int iterationStart= 0;
    double largestMNorm = 1e-50;
    double PERMITTIVITY_IRON = 2.22; // cobalt = 1.72, iron = 2.22, nickel = 0.6;
    double totalTime = 0;

    // #####################################        Protected Data Structures        ####################################

    // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins
    std::vector<double> exchangeVec;

    // Sites to be printed if printFixedSites is TRUE.
    std::list <int> fixedOutputSites;

    // Description missing
    std::vector<int> layerSpinPairs;

    // Description missing
    std::vector<int>    layerSpinsInChain;

    // Description missing
    std::vector<int>    layerTotalSpins;

    // Vectors containing magnetic components (m), along each axis, at the initial conditions for all spins. Leave as zero!

    std::vector<std::vector<std::vector<double>>> m0Nest;
    std::vector<std::vector<std::vector<double>>> m1Nest;
    std::vector<std::vector<std::vector<double>>> m2Nest;
};


#endif //SPINCHAINS_SYSTEMDATACONTAINER_H
