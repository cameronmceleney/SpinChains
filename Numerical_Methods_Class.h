#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include "CommonLibs.h"

class Numerical_Methods_Class {

private:
//  Dtype               Member Name                             //Comment

    double              _biasField = 0.1;                    // Bias field (T)
    double              _biasFieldDriving = 0;               // Driving field amplitude [T] (caution: papers often give in mT)
    double              _biasFieldDrivingScale = 2;             // The factor by which the driving field amplitude will be modulated
    double              _biasFieldDrivingUse;                   // Value to be used after bool statement. Will be either _biasFieldDrivingInit or _biasFieldDrivingShock
    double              _biasFieldDrivingInit;                  // Driving field amplitude [T] for use prior to shockwave. Commonly will be set to equal _biasFieldDriving
    double              _biasFieldDrivingShock;                 // Driving field amplitude [T] for the shockwave. Must be different to _biasFieldDriving to notice an effect
    std::vector<double> _chainJVals;                            // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins

    double              _drivingAngFreq;                        // Angular frequency of oscillatory driving field[rad*s^{-1}]
    double              _drivingFreq = 42.5 * 1e9;        // Frequency of oscillatory driving field [GHz] (f_d in literature) (default: 10 * 6.045 * 1e9)
    int                 _drivingRegionLHS;                      // The position of the spin which is leftmost in the driving region
    int                 _drivingRegionRHS;                      // The position of the spin which is rightmost in the driving region
    int                 _drivingRegionWidth;                    // Driving region width

    double              _gilbertConst = 1e-4;                   // Gilbert Damping Factor
    double              _gyroMagConst = 29.2E9 * 2 * M_PI;      // Gyromagnetic ratio (GHz/T). 29.2E9 is the numerical value of the gyromagetic ratio of the electron divided by 2pi
    bool                _hasShockWaveBegan;
    double              _magSat = 1.0;                          // Saturation Magnetisation (T). Note: 1A/m = 1.254uT. Must be in Telsa,
    double              _maxSimTime;                            // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time

    std::vector<double> _mxStartVal{0};                         // Vector containing magnetic components (m) along the x-axis (x) at the initial conditions for all spins
    std::vector<double> _myStartVal{0};                         // Vector containing magnetic components (m) along the y-axis (y) at the initial conditions for all spins
    std::vector<double> _mzStartVal{0};                         // Vector containing magnetic components (m) along the z-axis (z) at the initial conditions for all spins
    double              _mxInit = 0.001;                            // The initial value of the magnetic moment (m) along the x-direction. Normalised to mZInit (1.0 = 100% of mzInit)
    double              _myInit = 0;                            // The initial value of the magnetic moment (m) along the y-direction. Normalised to mZInit (1.0 = 100% of mzInit)
    double              _mzInit;                                // The initial value of the magnetic moment (m) along the z-direction.

    int                 _numberOfDataPoints = 100;             // How many data-points will be saved. Higher number gives greater precision, but drastically increases filesize. Default is 100.
    int                 _numberOfSpinPairs;                     // Number of pairs of spins in the chain. Used for array lengths and tidying notation
    bool                _shouldDebug = false;                   // Internal flag to indicate if debugging and output flags should be used, regardless of CMAKE build options
    double               _startIterVal = 0;                      // The iteration step that the program will begin at. Often set as zero
    double              _stepsize;                              // Stepsize between values
    double              _stepsizeHalf;                          // Separately defined to avoid repeated unnecessary calculations inside loops

    std::string         _stepsizeString;                        // Object to string conversation for value
    std::string         _stopIterString;                        // Object to string conversion for value
    double              _stopIterVal;                           // The maximum iteration step that the program will calculate to
    double              _totalTime = 0;                         // Analogous to a stopwatch in a physical experiment. This tracks for how long the experiment in the model has been simulated

public:
//  Dtype               Member Name                             // Comment
    void                NMSetup();
    void                RK2();
    void                RK2LLG();                               // Testing function to add nonlinearity test to original RK2 code
    void                RK2Shockwaves();                        // Testing function for shockwaves work
    void                StreamToString();
};

#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H