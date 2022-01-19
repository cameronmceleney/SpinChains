#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include "CommonLibs.h"

class Numerical_Methods_Class {

private:

    float               _biasField = 0.1;                               // Bias field (T)
    double              _biasFieldDrivingAmplitude = 3e-3;              // Driving field amplitude [T] (caution: papers often give in mT)
    std::vector<double> _chainExchangeValues;                           // Holds a linearly spaced array of values which describe the exchange interaction between neighbouring spins

    double              _drivingAngularFrequency;                       // Angular frequency [rad*s^{-1}]
    double              _drivingFrequency = 62.8;                       // Driving frequency [GHz] (f_d in literature)
    int                 _drivingRegionLHS;                              // The position of the spin which is leftmost in the driving region
    int                 _drivingRegionRHS;                              // The position of the spin which is rightmost in the driving region
    int                 _drivingRegionWidth;                            // Driving region width

    float               _gyroscopicMagneticConstant = 29.2E9 * 2 * M_PI;// Gyromagnetic ratio (GHz/T)
    double              _initialMagMomentX = 0;                         // The initial value of the magnetic moment (MM) along the x-direction
    double              _initialMagMomentY = 0;                         // The initial value of the magnetic moment (MM) along the y-direction
    double              _initialMagMomentZ;                             // The initial value of the magnetic moment (MM) along the z-direction

    float               _magnetisationSaturation = 1.0;                 // Saturation Magnetisation (T). Note: 1A/m = 1.254uT. Must be in Telsa,
    double              _maxSimulatedTime;                              // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time
    long                _startIterationValue = 0;                       // The minimum iteration step that the program will calculate to

    double              _stepSize;                                      // Stepsize between values
    double              _stepsizeHalf;                                  // Separately defined to avoid repeated unnecessary calculations inside loops

    std::string         _stepsizeString;                                // Object to string conversation for value
    std::string         _stopIterString;                                // Object to string conversion for value
    double              _stopIterationValue;                            // The maximum iteration step that the program will calculate to
    bool                _shouldDebug = false;                           // Internal flag to indicate if debugging and output flags should be used, regardless of CMAKE build options
    double              _totalTime = 0;                                 // Analogous to a stopwatch in a physical experiment. This tracks for how long the experiment in the model has been simulated



public:

    void                NMSetup();
    void                RK2();
    void                StreamToString();
};


#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H
