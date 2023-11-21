//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/MagnetisationDynamics.h"

MagnetisationDynamics::MagnetisationDynamics(SimulationParameters* sharedSimParams, 
                                             SimulationStates* sharedSimStates, 
                                             SimulationFlags* sharedSimFlags)
                                   
   : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}


double MagnetisationDynamics::MagneticMomentX(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mx;

    if ( _simFlags->shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mx = _simParams->gyroMagConst * (- (_simStates->gilbertVector[site] * hyMID * mxMID * myMID) + hyMID * mzMID - hzMID * (myMID
                + _simStates->gilbertVector[site] * mxMID * mzMID) + _simStates->gilbertVector[site] * hxMID * (pow(myMID,2) + pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mx = -1.0 * _simParams->gyroMagConst * (myMID * hzMID - mzMID * hyMID);
    }

    return mx;
}
double MagnetisationDynamics::MagneticMomentY(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double my;

    if ( _simFlags->shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        my = _simParams->gyroMagConst * (-(hxMID * mzMID) + hzMID * (mxMID - _simStates->gilbertVector[site] * myMID * mzMID) + _simStates->gilbertVector[site] * (hyMID * pow(mxMID,2) - hxMID * mxMID * myMID + hyMID * pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        my = _simParams->gyroMagConst * (mxMID * hzMID - mzMID * hxMID);
    }

    return my;
}
double MagnetisationDynamics::MagneticMomentZ(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mz;

    if ( _simFlags->shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mz = _simParams->gyroMagConst * (hxMID * myMID + _simStates->gilbertVector[site] * hzMID * (pow(mxMID,2) + pow(myMID,2)) - _simStates->gilbertVector[site]*hxMID*mxMID*mzMID - hyMID * (mxMID + _simStates->gilbertVector[site] * myMID * mzMID));
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
std::vector<double> MagnetisationDynamics::StochasticTerm(const int& site, const double &timeStep) {
    // Function to compute the stochastic term

    // Compute the standard deviation for the Gaussian noise
    double stddev = std::sqrt(2.0 * _simStates->gilbertVector[site] * SimulationParameters::BOLTZMANN_CONSTANT * _simParams->ambientTemperature / (_simParams->gyroMagConst * _simParams->satMag * timeStep));

    // Generate Gaussian noise for each direction
    double xi_x = GenerateGaussianNoise(0.0, stddev);
    double xi_y = GenerateGaussianNoise(0.0, stddev);
    double xi_z = GenerateGaussianNoise(0.0, stddev);

    return {xi_x, xi_y, xi_z};
}
std::vector<double> MagnetisationDynamics::ComputeStochasticTerm(const int& site, const double &timeStep) {
    // Function to compute the stochastic term
    std::vector<double> noise = StochasticTerm(site, timeStep);
    std::vector<double> stochasticField = {noise[0], noise[1], noise[2]};
    return stochasticField;
}
