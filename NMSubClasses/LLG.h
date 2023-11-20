//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_LLG_H
#define SPINCHAINS_LLG_H

#include "NMMethods.h"
#include "ILLG.h"

class MagnetisationDynamics : public NMMethods, public ILLG {
public:
    MagnetisationDynamics(std::shared_ptr<SystemDataContainer> data) : NMMethods(data) {}
    ~MagnetisationDynamics() = default;
public:
            /*
     * ################################################################################################################
     * #############    Terms to calculate the magnetic moments of the atoms (doesn't yet include sLLG)    ############
     * ################################################################################################################
     */

    // Simple version of magnetic moment x-component for single-layered systems
    double              MagneticMomentX (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) override;

    // Simple version of magnetic moment y-component for single-layered systems
    double              MagneticMomentY (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) override;

    // Simple version of magnetic moment z-component for single-layered systems
    double              MagneticMomentZ (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) override;

    // Full version of magnetic moment x-component function for multi-layered systems
    double              MagneticMomentX (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) override;

    // Full version of magnetic moment y-component function for multi-layered systems
    double              MagneticMomentY (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) override;

    // Full version of magnetic moment z-component function for multi-layered systems
    double              MagneticMomentZ (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID) override;

    /*
     * ################################################################################################################
     * ################    Terms to calculate the stochastic temperature contribution to sLLG (W.I.P)    ##############
     * ################################################################################################################
     */

    // Description missing
    double GenerateGaussianNoise(const double &mean, const double &stddev) override;

    // Description missing
    std::vector<double> StochasticTerm(const int& site, const double &timeStep) override;

    // Description missing
    std::vector<double> ComputeStochasticTerm(const int& site, const double &timeStep) override;
};


#endif //SPINCHAINS_LLG_H
