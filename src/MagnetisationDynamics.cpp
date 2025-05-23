//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/MagnetisationDynamics.h"

MagnetisationDynamics::MagnetisationDynamics(SimulationParameters* sharedSimParams, 
                                             SimulationStates* sharedSimStates, 
                                             SimulationFlags* sharedSimFlags)
                                   
   : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}

double MagnetisationDynamics::_checkIfDampingMapExists(const int& site) {
    // Function to check if the damping map exists

    if ( !_simFlags->hasGradientRegionForDamping ) {
        return _simStates->gilbertVector[site];
    }
    else {
        // Has a gradient map requested
        auto it = _simStates->dampingGradientMap.find(site);
        if (it != _simStates->dampingGradientMap.end()) { return it->second;}  // Found the site in the map
        else { return _simStates->gilbertVector[site]; }  // Site not found in the map
    }
}

double
MagnetisationDynamics::MagneticMomentX( const int &site, const double &mxMID, const double &myMID, const double &mzMID,
                                        const double &hxMID, const double &hyMID, const double &hzMID,
                                        const double &gilbertFactor ) {

    double mx;

    if ( _simFlags->shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mx = _simParams->gyroMagConst * (- (gilbertFactor * hyMID * mxMID * myMID) + hyMID * mzMID - hzMID * (myMID
                + gilbertFactor * mxMID * mzMID) + gilbertFactor * hxMID * (pow(myMID,2) + pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mx = -1.0 * _simParams->gyroMagConst * (myMID * hzMID - mzMID * hyMID);
    }

    return mx;
}
double
MagnetisationDynamics::MagneticMomentY( const int &site, const double &mxMID, const double &myMID, const double &mzMID,
                                        const double &hxMID, const double &hyMID, const double &hzMID,
                                        const double &gilbertFactor ) {

    double my;

    if ( _simFlags->shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        my = _simParams->gyroMagConst * (-(hxMID * mzMID) + hzMID * (mxMID - gilbertFactor * myMID * mzMID) + gilbertFactor * (hyMID * pow(mxMID,2) - hxMID * mxMID * myMID + hyMID * pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        my = _simParams->gyroMagConst * (mxMID * hzMID - mzMID * hxMID);
    }

    return my;
}
double
MagnetisationDynamics::MagneticMomentZ( const int &site, const double &mxMID, const double &myMID, const double &mzMID,
                                        const double &hxMID, const double &hyMID, const double &hzMID,
                                        const double &gilbertFactor ) {

    double mz;

    if ( _simFlags->shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mz = _simParams->gyroMagConst * (hxMID * myMID + gilbertFactor * hzMID * (pow(mxMID,2) + pow(myMID,2)) - gilbertFactor*hxMID*mxMID*mzMID - hyMID * (mxMID + gilbertFactor * myMID * mzMID));

    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mz = -1.0 * _simParams->gyroMagConst * (mxMID * hyMID - myMID * hxMID);
    }

    return mz;
}

double MagnetisationDynamics::MagneticMomentX(const int& site, const int& layer, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mx;

    if ( _simFlags->shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mx = _simParams->gyroMagConst * (- (_simStates->gilbertVectorMulti[layer][site] * hyMID * mxMID * myMID) + hyMID * mzMID - hzMID * (myMID + _simStates->gilbertVectorMulti[layer][site] * mxMID * mzMID) + _simStates->gilbertVectorMulti[layer][site] * hxMID * (pow(myMID,2) + pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mx = -1.0 * _simParams->gyroMagConst * (myMID * hzMID - mzMID * hyMID);
    }

    return mx;
}
double MagnetisationDynamics::MagneticMomentY(const int& site, const int& layer, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double my;

    if (_simFlags->shouldUseLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        my = _simParams->gyroMagConst * (-(hxMID * mzMID) + hzMID * (mxMID - _simStates->gilbertVectorMulti[layer][site] * myMID * mzMID) + _simStates->gilbertVectorMulti[layer][site] * (hyMID * pow(mxMID,2) - hxMID * mxMID * myMID + hyMID * pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        my = _simParams->gyroMagConst * (mxMID * hzMID - mzMID * hxMID);
    }

    return my;
}
double MagnetisationDynamics::MagneticMomentZ(const int& site, const int& layer, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mz;

    if (_simFlags->shouldUseLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mz = _simParams->gyroMagConst * (hxMID * myMID + _simStates->gilbertVectorMulti[layer][site] * hzMID * (pow(mxMID,2) + pow(myMID,2)) - _simStates->gilbertVectorMulti[layer][site]*hxMID*mxMID*mzMID - hyMID * (mxMID + _simStates->gilbertVectorMulti[layer][site] * myMID * mzMID));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mz = -1.0 * _simParams->gyroMagConst * (mxMID * hyMID - myMID * hxMID);
    }

    return mz;
}

double MagnetisationDynamics::GenerateGaussianNoise(const double &mean, const double &stddev) {
    // Function to generate random numbers from a Gaussian distribution
    static std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> distribution(mean, stddev);
    return distribution(generator);
}
std::vector<double>
MagnetisationDynamics::StochasticTerm( const int &site, const double &timeStep, const double &gilbertFactor ) {
    // Function to compute the stochastic term

    // Compute the standard deviation for the Gaussian noise
    double stddev = std::sqrt(2.0 * gilbertFactor * SimulationParameters::BOLTZMANN_CONSTANT * _simParams->ambientTemperature / (_simParams->gyroMagConst * _simParams->satMag * timeStep));

    // Generate Gaussian noise for each direction
    double xi_x = GenerateGaussianNoise(0.0, stddev);
    double xi_y = GenerateGaussianNoise(0.0, stddev);
    double xi_z = GenerateGaussianNoise(0.0, stddev);

    return {xi_x, xi_y, xi_z};
}
std::vector<double>
MagnetisationDynamics::ComputeStochasticTerm( const int &site, const double &timeStep, const double &gilbertFactor ) {
    // Function to compute the stochastic term
    std::vector<double> noise = StochasticTerm(site, timeStep, gilbertFactor);
    std::vector<double> stochasticField = {noise[0], noise[1], noise[2]};
    return stochasticField;
}
