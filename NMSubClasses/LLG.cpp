//
// Created by Cameron McEleney on 31/10/2023.
//

#include "LLG.h"

double MagnetisationDynamics::MagneticMomentX(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mxK;

    if (systemData->useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mxK = systemData->gyroMagConst * (- (systemData->gilbertVector[site] * hyMID * mxMID * myMID) + hyMID * mzMID - hzMID * (myMID
                + systemData->gilbertVector[site] * mxMID * mzMID) + systemData->gilbertVector[site] * hxMID * (pow(myMID,2) + pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mxK = -1.0 * systemData->gyroMagConst * (myMID * hzMID - mzMID * hyMID);
    }

    return mxK;
}
double MagnetisationDynamics::MagneticMomentY(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double myK;

    if (systemData->useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        myK = systemData->gyroMagConst * (-(hxMID * mzMID) + hzMID * (mxMID - systemData->gilbertVector[site] * myMID * mzMID) + systemData->gilbertVector[site] * (hyMID * pow(mxMID,2) - hxMID * mxMID * myMID + hyMID * pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        myK = systemData->gyroMagConst * (mxMID * hzMID - mzMID * hxMID);
    }

    return myK;
}
double MagnetisationDynamics::MagneticMomentZ(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mzK;

    if (systemData->useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mzK = systemData->gyroMagConst * (hxMID * myMID + systemData->gilbertVector[site] * hzMID * (pow(mxMID,2) + pow(myMID,2)) - systemData->gilbertVector[site]*hxMID*mxMID*mzMID - hyMID * (mxMID + systemData->gilbertVector[site] * myMID * mzMID));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mzK = -1.0 * systemData->gyroMagConst * (mxMID * hyMID - myMID * hxMID);
    }

    return mzK;
}

double MagnetisationDynamics::MagneticMomentX(const int& site, const int& layer, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mxK;

    if (systemData->useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mxK = systemData->gyroMagConst * (- (systemData->gilbertVectorMulti[layer][site] * hyMID * mxMID * myMID) + hyMID * mzMID - hzMID * (myMID + systemData->gilbertVectorMulti[layer][site] * mxMID * mzMID) + systemData->gilbertVectorMulti[layer][site] * hxMID * (pow(myMID,2) + pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mxK = -1.0 * systemData->gyroMagConst * (myMID * hzMID - mzMID * hyMID);
    }

    return mxK;
}
double MagnetisationDynamics::MagneticMomentY(const int& site, const int& layer, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double myK;

    if (systemData->useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        myK = systemData->gyroMagConst * (-(hxMID * mzMID) + hzMID * (mxMID - systemData->gilbertVectorMulti[layer][site] * myMID * mzMID) + systemData->gilbertVectorMulti[layer][site] * (hyMID * pow(mxMID,2) - hxMID * mxMID * myMID + hyMID * pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        myK = systemData->gyroMagConst * (mxMID * hzMID - mzMID * hxMID);
    }

    return myK;
}
double MagnetisationDynamics::MagneticMomentZ(const int& site, const int& layer, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mzK;

    if (systemData->useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mzK = systemData->gyroMagConst * (hxMID * myMID + systemData->gilbertVectorMulti[layer][site] * hzMID * (pow(mxMID,2) + pow(myMID,2)) - systemData->gilbertVectorMulti[layer][site]*hxMID*mxMID*mzMID - hyMID * (mxMID + systemData->gilbertVectorMulti[layer][site] * myMID * mzMID));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mzK = -1.0 * systemData->gyroMagConst * (mxMID * hyMID - myMID * hxMID);
    }

    return mzK;
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
    double stddev = std::sqrt(2.0 * systemData->gilbertVector[site] * systemData->BOLTZMANN_CONSTANT * systemData->ambientTemperature / (systemData->gyroMagConst * systemData->satMag * timeStep));

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
