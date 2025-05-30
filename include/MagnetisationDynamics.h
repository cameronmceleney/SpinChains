//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_MAGNETISATIONDYNAMICS_H
#define SPINCHAINS_MAGNETISATIONDYNAMICS_H

// C++ Standard Libraries
#include <cmath>
#include <random>

// C++ User Libraries (Containers)
#include "SimulationFlags.h"
#include "SimulationParameters.h"
#include "SimulationStates.h"

class MagnetisationDynamics {
private:
    SimulationParameters* _simParams; // Non-owning pointer to SimulationParameters
    SimulationStates* _simStates;
    SimulationFlags* _simFlags;
public:
    explicit MagnetisationDynamics(SimulationParameters* sharedSimParams,
                                   SimulationStates* sharedSimStates,
                                   SimulationFlags* sharedSimFlags);
    ~MagnetisationDynamics() = default;
public:
            /*
     * ################################################################################################################
     * #############    Terms to calculate the magnetic moments of the atoms (doesn't yet include sLLG)    ############
     * ################################################################################################################
     */

    // Simple version of magnetic moment x-component for single-layered systems
    double MagneticMomentX( const int &spin, const double &mxMID, const double &myMID, const double &mzMID,
                            const double &hxMID, const double &hyMID, const double &hzMID,
                            const double &gilbertFactor );

    // Simple version of magnetic moment y-component for single-layered systems
    double MagneticMomentY( const int &spin, const double &mxMID, const double &myMID, const double &mzMID,
                            const double &hxMID, const double &hyMID, const double &hzMID,
                            const double &gilbertFactor );

    // Simple version of magnetic moment z-component for single-layered systems
    double MagneticMomentZ( const int &spin, const double &mxMID, const double &myMID, const double &mzMID,
                            const double &hxMID, const double &hyMID, const double &hzMID,
                            const double &gilbertFactor );

    // Full version of magnetic moment x-component function for multi-layered systems
    double              MagneticMomentX (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    // Full version of magnetic moment y-component function for multi-layered systems
    double              MagneticMomentY (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    // Full version of magnetic moment z-component function for multi-layered systems
    double              MagneticMomentZ (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    /*
     * ################################################################################################################
     * ################    Terms to calculate the stochastic temperature contribution to sLLG (W.I.P)    ##############
     * ################################################################################################################
     */

    // Description missing
    double GenerateGaussianNoise(const double &mean, const double &stddev);

    // Description missing
    std::vector<double> StochasticTerm( const int &site, const double &timeStep, const double &gilbertFactor );

    // Description missing
    std::vector<double> ComputeStochasticTerm( const int &site, const double &timeStep, const double &gilbertFactor );

    double _checkIfDampingMapExists(const int& site);
};


#endif //SPINCHAINS_MAGNETISATIONDYNAMICS_H
