#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include "CommonLibs.h"

class Numerical_Methods_Class {

private:
//  Dtype               Member Name                             //Comment

    float               _biasField = 0.1;                       // Bias field (T)
    double              _biasFieldDriving = 3e-3;               // Driving field amplitude [T] (caution: papers often give in mT)
    std::vector<double> _chainJVals;                            // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins

    double              _drivingAngFreq;                        // Angular frequency of oscillatory driving field[rad*s^{-1}]
    double              _drivingFreq = 42.5e9;                  // Frequency of oscillatory driving field [GHz] (f_d in literature)
    int                 _drivingRegionLHS;                      // The position of the spin which is leftmost in the driving region
    double              _drivingRegionRHS;                      // The position of the spin which is rightmost in the driving region
    double              _drivingRegionWidth;                    // Driving region width

    float               _gyroMagConst = 29.2E9 * 2 * M_PI;      // Gyromagnetic ratio (GHz/T). 29.2E9 is the numerical value of the gyromagetic ratio of the electron divided by 2pi
    float               _magSat = 1.0;                          // Saturation Magnetisation (T). Note: 1A/m = 1.254uT. Must be in Telsa,
    double              _maxSimTime;                            // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time

    std::vector<double> _mxStartVal{0};                         // Vector containing magnetic components (m) along the x-axis (x) at the initial conditions for all spins
    std::vector<double> _myStartVal{0};                         // Vector containing magnetic components (m) along the y-axis (y) at the initial conditions for all spins
    std::vector<double> _mzStartVal{0};                         // Vector containing magnetic components (m) along the z-axis (z) at the initial conditions for all spins
    double              _mxInit = 0;                            // The initial value of the magnetic moment (m) along the x-direction
    double              _myInit = 0;                            // The initial value of the magnetic moment (m) along the y-direction
    double              _mzInit;                                // The initial value of the magnetic moment (m) along the z-direction

    int                 _numberOfSpinPairs;                      // Number of pairs of spins in the chain. Used for array lengths and tidying notation
    bool                _shouldDebug = false;                   // Internal flag to indicate if debugging and output flags should be used, regardless of CMAKE build options
    long                _startIterVal = 0;                      // The iteration step that the program will begin at. Often set as zero
    double              _stepsize;                              // Stepsize between values
    double              _stepsizeHalf;                          // Separately defined to avoid repeated unnecessary calculations inside loops

    std::string         _stepsizeString;                        // Object to string conversation for value
    std::string         _stopIterString;                        // Object to string conversion for value
    double              _stopIterVal;                           // The maximum iteration step that the program will calculate to
    double              _totalTime = 0;                         // Analogous to a stopwatch in a physical experiment. This tracks for how long the experiment in the model has been simulated

public:
//  Dtype               Member Name                                     //Comment
    void                NMSetup();
    void                RK2();
    void                StreamToString();
};

#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H