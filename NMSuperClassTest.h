//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMSUPERCLASSTEST_H
#define SPINCHAINS_NMSUPERCLASSTEST_H

// C++ Standard Library
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <memory>

// C++ Third Party Libraries
extern "C" {
    #include <fftw3.h>
}

// C++ User Libraries (General)
#include "linspace.h"
#include "SpinChainEigenSolverClass.h"
#include "CommonLibs.h"
#include "progressbar.hpp"

class NMSuperClassTest{
private:
    // ####################################            Private Methods            ###################################
    // ####################################            Protected Instances            ###################################


protected:
     // ####################################            Protected Instances            ###################################

//  Dtype               Member Name                                Variable docstring
     // ####################################            Protected Variables            ###################################
    static double               exchangeEnergyMin;                      // Minimum value of the exchange interaction energy (J) from the Heisenberg Hamiltonian, in [T].
    static double               exchangeEnergyMax;                      // Maximum value of the exchange interaction energy (J) from the Heisenberg Hamiltonian, in [T].

    static double              ambientTemperature;
    static double              anisotropyField;
    constexpr static double    BOLTZMANN_CONSTANT = 1.380649e-23;         // Boltzmann Constant [m^{2} kg s^{-2} K^{-1}].

    static double              drivingAngFreq;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
    static double              drivingFreq;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
    static int                 drivingRegionLhs;                         // The position of the spin which is leftmost in the driving region.
    static int                 drivingRegionRhs;                         // The position of the spin which is rightmost in the driving region.

    static int                 drivingRegionWidth;                       // Driving region width.
    static double              dynamicBiasField;                         // Driving field amplitude [T] (caution: papers often give in [mT]).
    static int                 forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
    static double              dipoleConstant;                           // Scaling factor which is constant across dipolar interaction calculations.

    static double              gilbertDamping;                           // Gilbert damping factor for main chain.
    static double              gilbertABCInner;                          // The lower Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the chain meets the ABC.
    static double              gilbertABCOuter;                          // The upper Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the ABC meets the pinned sites.
    static double              gyroMagConst;                             // Gyromagnetic ratio of an electron [GHz/T].

    static int                 iterationEnd;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    static int                 iterationStart;                       // The iteration step that the program will begin at. (Default: 0.0)
    static double              iterStartShock;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the maxSimTime.

    static double              iterEndShock;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the maxSimTime.
    static double              largestMNorm;                     // Computes sqrt(mxInit**2 + myInit**2 + mzInit**2). Initialised to be arbitrarily small.
    static double              maxSimTime;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [mxInit + myInit + mzInit]  CANNOT sum to greater than 1.0
    static double              PERMITTIVITY_IRON; // Magnetic moment in µB for iron

    static int                 numberNeighbours;
    static int                 numberOfDataPoints;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    static int                 numberOfSpinPairs;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    static int                 numSpinsDamped;                           // Number of spins in the damped regions (previously called _numGilbert).

    static int                 numSpinsInChain;                          // The number of spin sites in the spin chain to be simulated.
    static int                 systemTotalSpins;                         // The total number of spins in the system (chain plus ABCs).
    constexpr static double    PERM_FREESPACE = 1.25663706212e-6;         // Permeability of free space [H/m] (mu_0)
    static double              satMag;                                   // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)

    static double              shockwaveGradientTime;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    static double              shockwaveInitialStrength;                 // Initial strength of the shockwave before shockwaveScaling occurs. (Default: = dynamicBiasField)
    static double              shockwaveMax;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    static double              shockwaveScaling;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving

    static double              shockwaveStepsize;                        // Size of incremental increase in shockwave amplitude.
    static double              stepsize;                                 // Stepsize between values
    static double              stepsizeHalf;                             // Separately defined to avoid repeated unnecessary calculations inside loops
    static std::string         stepsizeString;                           // Object to string conversation for stepsize

    static std::string         stopIterString;                           // Object to string conversion for _stopIterVal
    static double              totalTime;                            // Analogous to a stopwatch in a physical experiment. This tracks the 'real time' of the simulation
    static int                 totalLayers;


    // #####################################        Protected Flags        ####################################
    static bool                centralDrive;                               // Drive from the centre of the chain if (true)
    static bool                driveAllLayers;
    static bool                dualDrive;                                  // Drive from both sides of the system
    static bool                hasShockwave;                               // Simulation contains a single driving bias field if (false).

    static bool                hasStaticDrive;                             // Selects (if true) whether drive has sinusoidal term
    static bool                isFm;
    static bool                isShockwaveOn;                              // Tests if the conditions to trigger a shockwave have been reached. Not to be altered by the user.
    static bool                isShockwaveAtMax;                           // Tests if the shockwave is at its maximum amplitude. Not to be altered by the user.

    static bool                lhsDrive;                                   // Drive from the LHS
    static bool                rhsDrive;                                   // Drive from the RHS
    static bool                printAllData;                               // Saves the m-component(s) of every spin at every iteration. WARNING: leads to huge output files.
    static bool                printFixedLines;                            // Saves m-component(s) of every spin at regular intervals. Total save points are set by numberOfDataPoints.
    static bool                printFixedSites;                            // Saves a discrete set of m-component(s) at regular intervals governed by numberOfDataPoints.

    static bool                shouldDriveCease;                           // Internal flag to indicate if the driving field should cut off at a given time.
    static bool                shouldTrackMValues;                         // Monitor the norm of all the m-values; if approx. 1.0 then the error is likely to be massive; discard that dataset.
    static bool                useLLG;                                     // Uses the Torque equation components if (false).
    static bool                useSLLG;

    static bool                useDipolar;
    static bool                useZeeman;
    static bool                useDemagIntense;                            // todo doesn't work
    static bool                useDemagFft;                                // todo doesn't work
    static bool                useMultilayer;
    static bool                debugFunc;

    // #####################################        Protected Data Structures        ####################################

    // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins
    static std::vector<double> exchangeVec;

    // Sites to be printed if printFixedSites is TRUE.
    static std::list <int> fixedOutputSites;

    // Description missing
    static std::vector<double> gilbertVector;

    // Description missing
    static std::vector<std::vector<double>> gilbertVectorMulti;

    // Computes sqrt(mxInit**2 + myInit**2 + mzInit**2). Initialised to be arbitrarily small.
    static std::vector<double> largestMNormMulti;

    // Description missing
    static std::vector<int> layerSpinPairs;

    // Description missing
    static std::vector<int>    layerSpinsInChain;

    // Description missing
    static std::vector<int>    layerTotalSpins;

    // Vectors containing magnetic components (m), along each axis, at the initial conditions for all spins. Leave as zero!

    // Magnetic moment along x-axis (m_x)
    static std::vector<double> mx0;

    // Magnetic moment along y-axis (m_y)
    static std::vector<double> my0;

    // Magnetic moment along z-axis (m_z)
    static std::vector<double> mz0;

    static std::vector<std::vector<std::vector<double>>> m0Nest;
    static std::vector<std::vector<std::vector<double>>> m1Nest;
    static std::vector<std::vector<std::vector<double>>> m2Nest;

public:
    NMSuperClassTest(); // Constructor
    ~NMSuperClassTest(); // Destructor (though not strictly necessary with smart pointers)

    virtual void childMethod() const;
    void callChildMethod() const;
    virtual ~Parent()
};


#endif //SPINCHAINS_NMSUPERCLASSTEST_H
