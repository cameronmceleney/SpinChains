//
// Created by Cameron Aidan McEleney on 19/11/2023.
//

#ifndef SPINCHAINS_SIMULATIONPARAMETERS_H
#define SPINCHAINS_SIMULATIONPARAMETERS_H

// C++ Standard Library

#define PERM_FREESPACE 1.25663706212e-6f  // Permeability of free space [H/m] (mu_0)
#define BOLTZMANN_CONSTANT 1.380649e-23f  // Boltzmann Constant [m^{2} kg s^{-2} K^{-1}].

struct SimulationParameters {
    // todo separate this out into several compositions (flags, data structures, etc)
    float PERMEABILITY_IRON = 2.22; // cobalt = 1.72, iron = 2.22, nickel = 0.6;

    // Sites to be printed if shouldPrintDiscreteSites is TRUE.

    float               gpu_current_time;
    float               ambientTemperature;
    float               anisotropyField;
    float              staticZeemanStrength;

    float               drivingAngFreq;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
    float               drivingFreq;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
    u_short              drivingRegionLhs;                         // The position of the spin which is leftmost in the driving region.
    u_short                 drivingRegionRhs;                         // The position of the spin which is rightmost in the driving region.

    u_int                 drivingRegionWidth;                       // Driving region width.
    float              oscillatingZeemanStrength;                         // Driving field amplitude [T] (caution: papers often give in [mT]).
    u_int                 forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
    float              dipoleConstant;                           // Scaling factor which is constant across dipolar interaction calculations.
    float              dmiConstant;
    float              exchangeStiffness;

    float              gilbertDamping;                           // Gilbert damping factor for main chain.
    float              gilbertABCInner;                          // The lower Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the chain meets the ABC.
    float              gilbertABCOuter;                          // The upper Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the ABC meets the pinned sites.
    float              gyroMagConst;                             // Gyromagnetic ratio of an electron [GHz/T].
    float              latticeConstant;

    float              spinPolarisation;                         // Spin polarisation of the spin current.
    float              spinTransferEfficiency;                   // Spin transfer efficiency of the spin current.

    u_int                 iterationEnd;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    u_int                 iterationStart = 0;                        // The iteration step that the program will begin at (Default: 0.0)
    float              iterStartShock;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the maxSimTime.

    float              iterEndShock;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the maxSimTime.
    float              largestMNorm = 1e-40;                     // Computes sqrt(_mx**2 + _my**2 + _mz**2) for each site at each moment to track any abnormalities. Initialises to be arbitrarily small
    float              maxSimTime;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [mxInit + myInit + mzInit]  CANNOT sum to greater than 1.0
    u_short                 numNeighbours;
    u_int                 numberOfDataPoints;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    u_int                 numberOfSpinPairs;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    u_int                 numSpinsInABC;                           // Number of spins in the damped regions (previously called _numGilbert).

    u_int                numSpinsInChain;                          // The number of spin sites in the spin chain to be simulated.
    u_int                 systemTotalSpins;                         // The total number of spins in the system (chain plus ABCs).
    float              satMag;                                   // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)

    float              shockwaveGradientTime;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    float              shockwaveInitialStrength;                 // Initial strength of the shockwave before shockwaveScaling occurs. (Default: = oscillatingZeemanStrength)
    float              shockwaveMax;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    float              shockwaveScaling;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving

    float              shockwaveStepsize;                        // Size of incremental increase in shockwave amplitude.
    float              stepsize;                                 // Stepsize between values
    float              stepsizeHalf;                             // Separately defined to avoid repeated unnecessary calculations inside loops
    u_short                 numLayers;
    float              recordingInterval;                        // Time between each data point recorded.
    u_short                 layerOfInterest;

    float exchangeEnergyMin;
    float exchangeEnergyMax;

    float totalTime = 0;

    u_short numSpinsDRPeak, numSpinsDRGradient;
    u_short dmiRegionLhs, dmiRegionRhs, numSpinsDmiPeak, numSpinsDmiGradient, numSpinsDmiWidth, dmiRegionOffset;
    u_short dampingRegionLhs, dampingRegionRhs, numSpinsDampingPeak, numSpinsDampingGradient, numSpinsDampingWidth, dampingRegionOffset;

    float dampingGradientPeak, dampingGradientGradient;

};

#endif //SPINCHAINS_SIMULATIONPARAMETERS_H
