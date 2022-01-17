#include "linspace.h"

void LinspaceClass::set_values(double intervalStart, double intervalEnd, int numberOfSamples, bool shouldIncludeEndpoint){
    // Value assignment is set separately to other class member functions in case value preprocessing of some sort is later required as added functionality

    // Start point of the interval.
    _intervalStart = intervalStart;

    // Stop point of the interval.
    _intervalEnd = intervalEnd;

    // Number (N) of evenly spaced samples to be returned over the given interval. Length of outputted linspace array is thus (N-1)
    _numberOfSamples = numberOfSamples;

    //If true, stop is appended to be the last sample. Otherwise, stop is not included.
    _shouldIncludeEndpoint = shouldIncludeEndpoint;
}

std::vector<double> LinspaceClass::generate_array() {

    if (_numberOfSamples == 0) {
        // if no range in inputted then empty vector is returned
        return _linspaceArray;
    }

    if (_numberOfSamples == 1) {
        // if range is one then start value is the only element of the array returned
        _linspaceArray.push_back(_intervalStart);
        return _linspaceArray;
    }

    // All other cases are dealt with here. Output will be of length 'num' if endpoint = True

    // 'delta' finds the stepsize between values
    double delta = (_intervalEnd - _intervalStart) / (_intervalStart - 1);

    for(int i=0; i < _numberOfSamples-1; ++i) {
        // Adds elements to list. Method is very similar to Euler's numerical method
        _linspaceArray.push_back(_intervalStart + delta * i);
    }

    if (_shouldIncludeEndpoint) {
        _linspaceArray.push_back(_intervalEnd);
    }

    return _linspaceArray;
}

std::vector<double> LinspaceClass::build_spinchain() {

    _spinchainArray.push_back(0); // Initialised with a zero to account for the (P-1)th spin
    _spinchainArray.insert(_spinchainArray.end(), _linspaceArray.begin(), _linspaceArray.end()); // Insert is faster for large values of numbers compared to push_back()
    _spinchainArray.push_back(0); // Appends a zero to the end to account for the exchange from the (N+1)th RHS spin

    return _spinchainArray;
}
