#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include <vector>
#include <string>

class Numerical_Methods_Class {

private:
    std::vector<double> _linspaceExchangeValues;

    std::string _fileName; //takes a custom filename input from the user

   double _PIVAL = 3.14159265359;

    int _numberOfSpins; // Number of sites (spins) in the chain
    int _spinStart = 0; // This is the site that the simulation will begin from
    int _drivingRegionLHSSpin = _spinStart; // The spin which is the leftmost in the driving region
    int _drivingRegionRHSSpin = _drivingRegionLHSSpin + 200; // The spin which is the rightmost in the driving region
    int _initalTime = 0; // Otherwise, known as t0, this is the start-time for the simulation
    int _startIterationValue = 0; //The minimum iteration step that the program will calculate to
    int _stopIterationValue; // The maximum iteration step that the program will calculate to

    double _stepsize; // Accepts float or scientific notation as input
    double _exchangeMinimum = 43.5;
    double _exchangeMaximum = 132;
    double _drivingFrequency = 67e9;
    double _drivingAngularFrequency = 2 * _PIVAL * _drivingFrequency; // angular freq (Hz)
    double _biasFieldDrivingAmplitude = 3e-3;
    double _maxSimulatedTime;     // notifies the user of the maximum simulation time. If the simulation time is too long, the user should simply force-exit the code

    float _biasField = 0.1; // bias field (T)
    float _magnetisationSaturation = 1.0; // Saturation Magnetisation (T). Note: 1A/m = 1.254uT. Must be in Telsa,
    float _gyroscopicMagneticConstant = 29.2E9 * 2 * _PIVAL; // gyromagnetic ratio (GHz/T)


public:

    void RK2(int numberSpins);
};


#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H
