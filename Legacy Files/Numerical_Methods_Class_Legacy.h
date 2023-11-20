#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include "../include/linspace.h"
#include "../include/SpinChainEigenSolverClass.h"
#include "../src/CommonLibs.h"
#include "../libs/progressbar.hpp"
#include <list>

class Numerical_Methods_Class {

private:
//  Dtype               Member Name                                Variable docstring

    double              _anisotropyField;
    double              _drivingAngFreq;                        // Angular frequency of oscillatory driving field [rad*s^{-1}].
    double              _drivingFreq;                           // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
    int                 _drivingRegionLHS;                      // The position of the spin which is leftmost in the driving region.
    int                 _drivingRegionRHS;                      // The position of the spin which is rightmost in the driving region.

    int                 _drivingRegionWidth;                    // Driving region width.
    double              _dynamicBiasField;                      // Driving field amplitude [T] (caution: papers often give in [mT]).
    std::list <int>     _fixed_output_sites;                    // Sites to be printed if _printFixedSites is TRUE.
    int                 _forceStopAtIteration;                  // Legacy breakpoint variable. Set as a -ve value to deactivate.
    double              _gilbertConst;                          // Gilbert Damping Factor.

    double              _gilbertLower;                          // The lower boundary for the damped regions at either end of the spinchain.
    double              _gilbertUpper;                          // The upper boundary for the damped regions at either end of the spinchain.
    double              _gyroMagConst;                          // Gyromagnetic ratio of an electron [GHz/T].
    int                 _iterationEnd;                          // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    int                 _iterationStart = 0;                    // The iteration step that the program will begin at. (Default: 0.0)

    double              _iterStartShock;                        // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the _maxSimTime.
    double              _iterEndShock;                          // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the _maxSimTime.
    double              _largestMNorm = 1e-50;                  // Computes sqrt(_mxInit**2 + _myInit**2 + _mzInit**2). Initialised to be arbitrarily small.
    double              _magSat = 1.0;                          // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)
    double              _maxSimTime;                            // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [_mxInit + _myInit + _mzInit]  CANNOT sum to greater than 1.0
    double              _mxInit = 0.0;                          // x-direction. (Default: 0.0)
    double              _myInit = 0.0;                          // y-direction. (Default: 0.0)
    double              _mzInit = _magSat;                      // z-direction. (Default: _magSat = 1.0)

    int                 _numberOfDataPoints;                    // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    int                 _numberOfSpinPairs;                     // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    int                 _numSpinsDamped;                        // Number of spins in the damped regions (previously called _numGilbert).
    int                 _numSpinsInChain;                       // The number of spin sites in the spin chain to be simulated.
    double              _recordingInterval;
    double              _regionScaling = 0.05;                  // Calculate _drivingRegionWidth as a fraction [0.0, 1.0] of _numSpinsInChain.

    double              _shockwaveGradientTime;                 // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    double              _shockwaveInitialStrength;              // Initial strength of the shockwave before _shockwaveScaling occurs. (Default: = _dynamicBiasField)
    double              _shockwaveMax;                          // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    double              _shockwaveScaling;                      // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving
    double              _shockwaveStepsize;                     // Size of incremental increase in shockwave amplitude.

    double              _stepsize;                              // Stepsize between values
    double              _stepsizeHalf;                          // Separately defined to avoid repeated unnecessary calculations inside loops
    std::string         _stepsizeString;                        // Object to string conversation for _stepsize
    std::string         _stopIterString;                        // Object to string conversion for _stopIterVal
    double              _totalTime = 0;                         // Analogous to a stopwatch in a physical experiment. This tracks the 'real' time' of the simulation

    // ######## Booleans and other tests ########
    bool                _centralDrive;                          // Drive from the centre of the chain if (true)
    bool                _dualDrive;                             // Drive from both sides of the system
    bool                _hasShockwave;                          // Simulation contains a single driving bias field if (false).
    bool                _hasStaticDrive;                        // Selects (if true) whether drive has sinusoidal term

    bool                _isFM;
    bool                _isShockwaveOn = false;                 // Tests if the conditions to trigger a shockwave have been reached. Not to be altered by the user.
    bool                _isShockwaveAtMax = false;              // Tests if the shockwave is at its maximum amplitude. Not to be altered by the user.
    bool                _lhsDrive;                              // Drive from the RHS if (false)

    bool                _printAllData;                          // Saves the m-component(s) of every spin at every iteration. WARNING: leads to huge output files.
    bool                _printFixedLines;                       // Saves m-component(s) of every spin at regular intervals. Total savepoints are set by _numberOfDataPoints.
    bool                _printFixedSites;                       // Saves a discrete set of m-component(s) at regular intervals governed by _numberOfDataPoints.
    bool                _shouldDriveCease;                      // Internal flag to indicate if the driving field should cut-off at a given time.

    bool                _shouldTrackMValues;                    // Monitor the norm of all the m-values; if approx. 1.0 then the error is likely to be massive; discard that dataset.
    bool                _useLLG;                                // Uses the Torque equation components if (false).

    // ######## Private Functions ########
    std::vector<double> _exchangeVec;                           // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins
    std::vector<double> _gilbertVector{0};

    // Vectors containing magnetic components (m), along each axis, at the initial conditions for all spins. Leave as zero!
    std::vector<double> _mx0{0};                                // x-axis (x)
    std::vector<double> _my0{0};                                // y-axis (y)
    std::vector<double> _mz0{0};                                // z-axis (z)

    // Private functions
    void                FinalChecks();
    void                SetShockwaveConditions();
    void                SetDampingRegion();
    void                SetDrivingRegion();
    void                SetExchangeVector();
    void                SetInitialMagneticMoments();

    void                CreateColumnHeaders(std::ofstream &outputFileName);
    void                CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata=false);
    void                CreateMetadata(bool print_end_time=false);
    void                InformUserOfCodeType(const std::string& nameNumericalMethod);
    void                PrintVector(std::vector<double> &vectorToPrint, bool shouldExitAfterPrint);
    void                SaveDataToFile(std::ofstream &outputFileName, std::vector<double> &arrayToWrite, int &iteration);
    void                TestShockwaveConditions(double iteration);

public:
//  Dtype               Member Name                                Variable docstring
    void                NMSetup();                              // Assignment of all values required for the simulation
    void                RK2OriginalFM();                        // Original RK2 method that I wrote. Left here as legacy working version of RK2. NOTE: heavily outdated
    void                RK2MidpointFM();                        // Evaluate the given system, using the Runge-Kutta (2nd Order) method, for a ferromagnetic material

    void                RK2MidpointFMForTesting();              // Debugging case for RK-2 (ferromagnetic) where a large number of output statements allow for all used parameters to be tracked throughout the simulation.
    void                RK4MidpointFM();                        // OUTDATED: Evaluate the given system, using the Runge-Kutta (4nd Order) method, for a ferromagnetic material
    void                RK2MidpointAFM();                       // Evaluate the given system, using the Runge-Kutta (2nd Order) method, for an anti-ferromagnetic material
};

#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H