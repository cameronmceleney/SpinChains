//
// Created by Cameron Aidan McEleney on 20/11/2023.
//

#ifndef SPINCHAINS_ILLG_H
#define SPINCHAINS_ILLG_H

#include <vector>

class ILLG {
public:
        // Simple version of magnetic moment x-component for single-layered systems
    virtual double              MagneticMomentX (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) = 0;

    // Simple version of magnetic moment y-component for single-layered systems
    virtual double              MagneticMomentY (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) = 0;

    // Simple version of magnetic moment z-component for single-layered systems
    virtual double              MagneticMomentZ (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) = 0;

    // Full version of magnetic moment x-component function for multi-layered systems
    virtual double              MagneticMomentX (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) = 0;

    // Full version of magnetic moment y-component function for multi-layered systems
    virtual double              MagneticMomentY (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) = 0;

    // Full version of magnetic moment z-component function for multi-layered systems
    virtual double              MagneticMomentZ (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) = 0;

    /*
     * ################################################################################################################
     * ################    Terms to calculate the stochastic temperature contribution to sLLG (W.I.P)    ##############
     * ################################################################################################################
     */

    // Description missing
    virtual double GenerateGaussianNoise(const double &mean, const double &stddev) = 0;

    // Description missing
    virtual std::vector<double> StochasticTerm(const int& site, const double &timeStep) = 0;

    // Description missing
    virtual std::vector<double> ComputeStochasticTerm(const int& site, const double &timeStep) = 0;

    ~ILLG() = default;
};
#endif //SPINCHAINS_ILLG_H
