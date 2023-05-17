#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include "linspace.h"
#include "SpinChainEigenSolverClass.h"
#include "CommonLibs.h"
#include "progressbar.hpp"
#include <chrono>
#include <iomanip>
#include <list>
#include <map>


class Numerical_Methods_Class {

private:
//  Dtype               Member Name                                Variable docstring

    double              _anisotropyField;
    double              _drivingAngFreq;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
    double              _drivingFreq;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
    int                 _drivingRegionLHS;                         // The position of the spin which is leftmost in the driving region.
    int                 _drivingRegionRHS;                         // The position of the spin which is rightmost in the driving region.

    int                 _drivingRegionWidth;                       // Driving region width.
    double              _dynamicBiasField;                         // Driving field amplitude [T] (caution: papers often give in [mT]).
    std::list <int>     _fixed_output_sites;                       // Sites to be printed if _printFixedSites is TRUE.
    int                 _forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
    double              _gilbertConst;                             // Gilbert Damping Factor.

    double              _gilbertLower;                             // The lower boundary for the damped regions at either end of the spinchain.
    double              _gilbertUpper;                             // The upper boundary for the damped regions at either end of the spinchain.
    double              _gyroMagConst;                             // Gyromagnetic ratio of an electron [GHz/T].
    int                 _iterationEnd;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    int                 _iterationStart = 0;                       // The iteration step that the program will begin at. (Default: 0.0)

    double              _iterStartShock;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the _maxSimTime.
    double              _iterEndShock;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the _maxSimTime.
    double              _largestMNorm = 1e-50;                     // Computes sqrt(_mxInit**2 + _myInit**2 + _mzInit**2). Initialised to be arbitrarily small.
    double              _magSat = 1.0;                             // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)
    double              _maxSimTime;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [_mxInit + _myInit + _mzInit]  CANNOT sum to greater than 1.0
    double              _mxInit = 0.0;                             // x-direction. (Default: 0.0)
    double              _myInit = 0.0;                             // y-direction. (Default: 0.0)
    double              _mzInit = _magSat;                         // z-direction. (Default: _magSat = 1.0)

    int                 _numberOfDataPoints;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    int                 _numberOfSpinPairs;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    int                 _numSpinsDamped;                           // Number of spins in the damped regions (previously called _numGilbert).
    int                 _numSpinsInChain;                          // The number of spin sites in the spin chain to be simulated.
    double              _recordingInterval;

    double              _permFreeSpace = 1.25663706212e-6;         // Permeability of free space [H/m] (mu_0)
    double              _shockwaveGradientTime;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    double              _shockwaveInitialStrength;                 // Initial strength of the shockwave before _shockwaveScaling occurs. (Default: = _dynamicBiasField)
    double              _shockwaveMax;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    double              _shockwaveScaling;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving
    double              _shockwaveStepsize;                        // Size of incremental increase in shockwave amplitude.

    double              _stepsize;                                 // Stepsize between values
    double              _stepsizeHalf;                             // Separately defined to avoid repeated unnecessary calculations inside loops
    std::string         _stepsizeString;                           // Object to string conversation for _stepsize
    std::string         _stopIterString;                           // Object to string conversion for _stopIterVal
    double              _totalTime = 0;                            // Analogous to a stopwatch in a physical experiment. This tracks the 'real time' of the simulation

    // ######## Booleans including tests ########
    bool                _centralDrive;                             // Drive from the centre of the chain if (true)
    bool                _dualDrive;                                // Drive from both sides of the system
    bool                _hasShockwave;                             // Simulation contains a single driving bias field if (false).
    bool                _hasStaticDrive;                           // Selects (if true) whether drive has sinusoidal term

    bool                _isFM;
    bool                _isShockwaveOn = false;                    // Tests if the conditions to trigger a shockwave have been reached. Not to be altered by the user.
    bool                _isShockwaveAtMax = false;                 // Tests if the shockwave is at its maximum amplitude. Not to be altered by the user.
    bool                _lhsDrive;                                 // Drive from the RHS if (false)

    bool                _printAllData;                             // Saves the m-component(s) of every spin at every iteration. WARNING: leads to huge output files.
    bool                _printFixedLines;                          // Saves m-component(s) of every spin at regular intervals. Total save points are set by _numberOfDataPoints.
    bool                _printFixedSites;                          // Saves a discrete set of m-component(s) at regular intervals governed by _numberOfDataPoints.
    bool                _shouldDriveCease;                         // Internal flag to indicate if the driving field should cut off at a given time.

    bool                _shouldTrackMValues;                       // Monitor the norm of all the m-values; if approx. 1.0 then the error is likely to be massive; discard that dataset.
    bool                _useLLG;                                   // Uses the Torque equation components if (false).
    bool                _useDipolar;
    bool                _useBilayer;

    // ######## Private Functions ########
    std::vector<double> _exchangeVec;                              // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins
    std::vector<double> _gilbertVector{0};

    // Vectors containing magnetic components (m), along each axis, at the initial conditions for all spins. Leave as zero!
    std::vector<double> _mx0{0};                                   // x-axis (x)
    std::vector<double> _my0{0};                                   // y-axis (y)
    std::vector<double> _mz0{0};                                   // z-axis (z)

    // Private functions
    void                FinalChecks();
    void                SetShockwaveConditions();
    void                SetDampingRegion();
    void                SetDrivingRegion();
    void                SetExchangeVector();
    void                SetInitialMagneticMoments();

    void                SetInitialMagneticMomentsMultilayer(std::vector<std::vector<std::vector<double>>>& nestedNestedVector,
                                                            int layer, double mxInit, double myInit, double mzInit);  // test
    std::vector<std::vector<std::vector<double>>> initializeNestedNestedVector(int numSites, bool includeEnd);
    void                SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration);

    void                CreateColumnHeaders(std::ofstream &outputFileName);
    void                CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata=false);
    void                CreateMetadata(bool print_end_time=false);
    void                InformUserOfCodeType(const std::string& nameNumericalMethod);
    void                PrintVector(std::vector<double> &vectorToPrint, bool shouldExitAfterPrint);
    void                PrintNestedNestedVector(std::vector<std::vector<std::vector<double>>> nestedNestedVector);
    void                SaveDataToFile(std::ofstream &outputFileName, std::vector<double> &arrayToWrite, int &iteration);
    void                TestShockwaveConditions(double iteration);

    std::vector<double> DipoleDipoleCoupling(std::vector<double> mxTerms, std::vector<double> myTerms,
                                             std::vector<double> mzTerms, std::vector<int> sitePositions);
    double              EffectiveFieldX (const int& site, const double& mxLHS, const double& mxMID,
                                         const double& mxRHS, const double& dipoleTerm, const double& current_time);
    double              EffectiveFieldY (const int& site, const double& myLHS, const double& myMID, const double& myRHS,
                                         const double& dipoleTerm);
    double              EffectiveFieldZ (const int& site, const double& mzLHS, const double& mzMID, const double& mzRHS,
                                         const double& dipoleTerm);

    double              MagneticMomentX (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);
    double              MagneticMomentY (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);
    double              MagneticMomentZ (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

public:
//  Dtype               Member Name                                Variable docstring
    void                NMSetup();                                 // Assignment of all values required for the simulation
    void                SolveRK2();                                // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2Bilayer();                         // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2BilayerTest();                                // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint
};

#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H
