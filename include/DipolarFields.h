//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_DIPOLARFIELDS_H
#define SPINCHAINS_DIPOLARFIELDS_H

// C++ Standard Library
#include <cstring>
#include <iostream>

// C++ Third Party Libraries
extern "C" {
    #include <fftw3.h>
}

// C++ User Libraries (General)
#include "CommonLibs.h"
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"


class DipolarFields {
private:
    SimulationParameters* _simParams; // Non-owning pointer to SimulationParameters
    SimulationStates* _simStates;
    SimulationFlags* _simFlags;
public:
    explicit DipolarFields(SimulationParameters* sharedSimParams,
                           SimulationStates* sharedSimStates,
                           SimulationFlags* sharedSimFlags);
    ~DipolarFields() = default;
public:
    // ####################################            Define Public Variables            ###################################

    // ####################################            Define Public Functions            ###################################
        // Description missing
    double dipolarKernel1D(const int& originSite, const int& influencingSite, const std::string& component);
    double dipolarKernel3D(const int& originSite, const int& influencingSite, const double& A, const double& alpha);
    void DipolarInteraction1D(std::vector<double> inMxTerms, std::vector<double>& outDipoleX);

    std::vector<double> DipolarInteractionClassic(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                     std::vector<double> mzTerms, std::vector<int> sitePositions);

    // Description missing
    std::vector<double> DipolarInteractionIntralayer(std::vector<std::vector<double>>& mTerms,
                                                     int& currentSite, const int& currentLayer=0,
                                                     const double& exchangeStiffness=5.3e-17);

    // Description missing
    std::vector<double> DipolarInteractionInterlayer(std::vector<std::vector<double>>& mTermsLayer1,
                                                     std::vector<std::vector<double>>& mTermsLayer2, int& currentSite,
                                                     const int& currentLayer, const int& otherLayer);

    // Description missing
    std::vector<double> DipolarInteractionInterlayerAll(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2,
                                                        int& currentSite, const int& currentLayer, const int& otherLayer,
                                                        double& exchangeStiffness, double& interlayerExchange);

    // Description missing
    std::vector<double> DipolarInteractionInterlayerAdjacent(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2, int& numNeighbours,
                                                        int& currentSite, const int& currentLayer,
                                                        double& exchangeStiffness, double& interlayerExchange);

    std::vector<double> DipolarInteractionIntralayerDebug(std::vector<std::vector<double>>& mTerms, int& numNeighbours,
                                              int& currentSite, const int& currentLayer = 0);
    // Description missing
    std::vector<double> DipolarInteractionInterlayerDebug(std::vector<std::vector<double>>& mTermsChain1,
                                                     std::vector<std::vector<double>>& mTermsChain2, int& numNeighbours,
                                                     int& currentSite, const int& currentLayer = 0);
};


#endif //SPINCHAINS_DIPOLARFIELDS_H
