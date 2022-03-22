#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include "linspace.h"
#include "SpinChainEigenSolverClass.h"
#include "CommonLibs.h"
#include "progressbar.hpp"

class Numerical_Methods_Class {

private:
//  Dtype               Member Name                             //Comment

    double              _biasField = 0.1;                       // Applied bias field [T]. Often written as H_0 or H_static in literature
    double              _biasFieldDriving;                      // Driving field amplitude [T] (caution: papers often give in mT)
    double              _shockwaveScaling;                      // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving

    std::vector<double> _chainJVals;                            // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins

    double              _drivingAngFreq;                        // Angular frequency of oscillatory driving field[rad*s^{-1}]
    double              _drivingFreq;                           // Frequency of oscillatory driving field [GHz] (f_d in literature) (default: 10 * 6.045 * 1e9)
    int                 _drivingRegionLHS;                      // The position of the spin which is leftmost in the driving region
    int                 _drivingRegionRHS;                      // The position of the spin which is rightmost in the driving region
    int                 _drivingRegionWidth;                    // Driving region width

    double              _gilbertConst = 1e-4;                   // Gilbert Damping Factor
    double              _gyroMagConst = 29.2E9 * 2 * M_PI;      // Gyromagnetic ratio of an electron [GHz/T].
    double              _iterToBeginShockwave;
    double              _linearFMR;
    double              _magSat = 1.0;                          // Saturation Magnetisation [T]. Note: 1A/m = 1.254uT. Must be in Telsa,
    double              _maxSimTime;                            // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time

    // Vectors containing magnetic components (m), along each axis, at the initial conditions for all spins
    std::vector<double> _mxStartVal{0};                         // x-axis (x)
    std::vector<double> _myStartVal{0};                         // y-axis (y)
    std::vector<double> _mzStartVal{0};                         // z-axis (z)

    /* The initial value of the magnetic moment (m) along each axis.
     * When setting [_mxInit, _myInit, _mzInit] note that they CANNOT sum to greater than 1.0
     */
    double              _mxInit = 0.0;                          // x-direction. Default is (0.0)
    double              _myInit = 0.0;                          // y-direction. Default is (0.0)
    double              _mzInit = _magSat;                      // z-direction. Default is (1.0)

    int                 _numberOfDataPoints;                    // How many data-points will be saved in the output file. Higher number gives greater precision, but drastically increases filesize. Default is 1000.
    int                 _numberOfSpinPairs;                     // Number of pairs of spins in the chain. Used for array lengths and tidying notation
    int                 _startIterVal = 0;                      // The iteration step that the program will begin at. Often set as zero
    double              _stepsize;                              // Stepsize between values
    double              _stepsizeHalf;                          // Separately defined to avoid repeated unnecessary calculations inside loops

    std::string         _stepsizeString;                        // Object to string conversation for value
    std::string         _stopIterString;                        // Object to string conversion for value
    int                 _stopIterVal;                           // The maximum iteration step that the program will calculate to
    double              _totalTime = 0;                         // Analogous to a stopwatch in a physical experiment. This tracks for how long the experiment in the model has been simulated

    bool                _useLLG;
    bool                _hasShockwave;
    bool                _isShockwaveAlreadyOn = false;
    bool                _shouldDebug = false;                   // Internal flag to indicate if debugging and output flags should be used, regardless of CMAKE build options
    bool                _saveAllSpins;
    bool                _onlyShowFinalState;

    // Private functions
    void                CreateColumnHeaders(std::ofstream &outputFileName, bool &areAllSpinBeingSaved, bool &onlyShowFinalState);
    void                CreateFileHeader(std::ofstream &outputFileName, bool &areAllSpinBeingSaved, bool &onlyShowFinalState);
    void                DebugOptions(std::vector<double> mxNextVal, std::vector<double> myNextVal, std::vector<double> mzNextVal, int spin, long iterationIndex);
    void                InformUserOfCodeType();
    void                SaveDataToFile(bool &areAllSpinBeingSaved, std::ofstream &outputFileName,
                                       std::vector<double> &arrayToWrite, int &iteration, bool &onlyShowFinalState);
    void                SetShockwaveConditions();
    void                SetupVectors();
    void                StreamToString();

public:
//  Dtype               Member Name                             // Comment
    void                NMSetup();
    void                RK2();
    void                RK2LLG();                               // Testing function to add nonlinearity test to original RK2 code
};

#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H