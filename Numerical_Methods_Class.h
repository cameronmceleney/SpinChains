#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include <vector>
#include <string>
#include <cmath>

class Numerical_Methods_Class {

private:
    std::vector<double> _chainExchangeValues; // Holds a linearly spaced array of values which describe the exchange interaction between neighbouring spins

    std::string _fileName; //takes a custom filename input from the user

   float _piconstantval = 3.14159265359;

    int _numberOfSpins; // Number of sites (spins) in the chain
    int _drivingRegionLHS; // The position of the spin which is leftmost in the driving region
    int _drivingRegionWidth = 200;
    int _drivingRegionRHS; // The position of the spin which is leftmost in the driving region
    long _startIterationValue = 0; //The minimum iteration step that the program will calculate to
    double _stopIterationValue; // The maximum iteration step that the program will calculate to

    double _stepsize; // Accepts float or scientific notation as input
    double _stepsizeHalf; // Separately defined to avoid repeated unnecessary calculations inside loops
    double _exchangeMinimum = 43.5;
    double _exchangeMaximum = 132;
    double _drivingFrequency = 42.5e9;
    double _drivingAngularFrequency = 2 * _piconstantval * _drivingFrequency; // angular freq (Hz)
    double _biasFieldDrivingAmplitude = 3e-3;
    double _maxSimulatedTime;     // notifies the user of the maximum simulation time. If the simulation time is too long, the user should simply force-exit the code
    float _magnetisationSaturation = 1.0; // Saturation Magnetisation (T). Note: 1A/m = 1.254uT. Must be in Telsa,
    double _initialMagMomentX = 0; // The initial value of the magnetic moment (MM) along the x-direction
    double _initialMagMomentY = 0; // The initial value of the magnetic moment (MM) along the y-direction
    double _initialMagMomentZ = _magnetisationSaturation; // The initial value of the magnetic moment (MM) along the z-direction
    double _totalTime = 0; // Analogous to a stopwatch in a physical experiment. This tracks for how long the experiment in the model has been simulated

    float _biasField = 0.1; // bias field (T)
    float _gyroscopicMagneticConstant = 29.2E9 * 2 * M_PI; // gyromagnetic ratio (GHz/T)

    bool _shouldDebug = false; // Internal flag to indicate if debugging and output flags should be used, regardless of CMAKE build options
    std::string _stepsizeString;
    std::string _stopIterString;

public:

    void RK2(int numberSpins);
    void StreamToString();
};


#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H
