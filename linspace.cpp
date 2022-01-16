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

std::vector<double> LinspaceClass::generate_array()
{
    std::vector<double> linspace_array;

    if (_numberOfSamples == 0) {
        // if no range in inputted then empty vector is returned
        return linspace_array;
    }

    if (_numberOfSamples == 1) {
        // if range is one then start value is the only element of the array returned
        linspace_array.push_back(_intervalStart);
        return linspace_array;
    }

    // All other cases are dealt with here. Output will be of length 'num' if endpoint = True

    // 'delta' finds the stepsize between values
    double delta = (_intervalEnd - _intervalStart) / (_intervalStart - 1);

    for(int i=0; i < _numberOfSamples-1; ++i) {
        // Adds elements to list. Method is very similar to Euler's numerical method
        linspace_array.push_back(_intervalStart + delta * i);
    }

    if (_shouldIncludeEndpoint) {
        linspace_array.push_back(_intervalEnd);
    }

    return linspace_array;
}

