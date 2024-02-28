#include "linspace.h"

// Setter function for linspace class
void
LinspaceClass::set_values( double intervalStart, double intervalEnd, int numberOfSamples, bool shouldIncludeEndpoint,
                           bool shouldBuildSpinChain, bool hasExclusiveBounds ) {

    _intervalStart = intervalStart;

    _intervalEnd = intervalEnd;

    _numberOfSamples = numberOfSamples;

    _shouldIncludeEndpoint = shouldIncludeEndpoint;

    _shouldBuildSpinChain = shouldBuildSpinChain;

    _hasExclusiveBounds = hasExclusiveBounds;

    if (_hasExclusiveBounds) {
        _shouldIncludeEndpoint = false; // To ensure that the range is over (_intervalStart, _intervalEnd) instead of (_intervalStart, _intervalEnd]
        _numberOfSamples += 2;
    }
}

// Getter function for linspace class
std::vector<double> LinspaceClass::generate_array() {
    _linspaceArray.clear(); // Required when using multi-layers as otherwise previous damping regions are retained

    if ( _numberOfSamples == 0 ) {
        // If no range in inputted then empty vector is returned
        if ( _shouldBuildSpinChain ) {
            _linspaceArray.push_back(0);
            _linspaceArray.push_back(0);
        }
        return _linspaceArray;
    }

    if ( _numberOfSamples == 1 ) {
        // If range is one then start value is the only element of the array returned
        if ( _shouldBuildSpinChain ) {
            _linspaceArray.push_back(0);
            _linspaceArray.push_back(_intervalStart);
            _linspaceArray.push_back(0);
        } else
            _linspaceArray.push_back(_intervalStart);

        return _linspaceArray;
    }

    // All other cases are dealt with here. Output will be of length 'num' if endpoint = True

    // 'delta' finds the stepsize between values
    double delta = (_intervalEnd - _intervalStart) / (_numberOfSamples - 1);
    int iStart = 0;
    if ( _hasExclusiveBounds ) { iStart = 1; }
    for ( int i = iStart; i < _numberOfSamples - 1; ++i ) {
        // Adds elements to list. Method is very similar to Euler's numerical method
        _linspaceArray.push_back(_intervalStart + delta * i);
    }

    if ( _shouldIncludeEndpoint ) { _linspaceArray.push_back(_intervalEnd); }

    if ( _shouldBuildSpinChain ) { build_spinchain(); }

    return _linspaceArray;
}

/* Special Getter function for linspace class. Must be used AFTER set_values() and generate_array() are called from
 * linspace class. This function appends zeros to start&end of array to represent spins at the end of the chain */
void LinspaceClass::build_spinchain() {

    // Initialised with a zero to account for the (P-1)th spin
    std::vector<double> spinChainArray = {0};

    // Insert is faster for large values of numbers compared to push_back()
    spinChainArray.insert(spinChainArray.end(), _linspaceArray.begin(), _linspaceArray.end());

    // Appends a zero to the end to account for the exchange from the (N+1)th RHS spin
    spinChainArray.push_back(0);

    _linspaceArray.clear();
    _linspaceArray = spinChainArray;
}

std::vector<double> LinspaceClass::build_spinchain_explicit( std::vector<double> linspaceArrayIn ) {

    // Initialised with a zero to account for the (P-1)th spin
    std::vector<double> spinChainArray = {0};

    // Insert is faster for large values of numbers compared to push_back()
    spinChainArray.insert(spinChainArray.end(), linspaceArrayIn.begin(), linspaceArrayIn.end());

    // Appends a zero to the end to account for the exchange from the (N+1)th RHS spin
    spinChainArray.push_back(0);

    return spinChainArray;
}
