//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_DIPOLARFIELD_H
#define SPINCHAINS_DIPOLARFIELD_H

#include "NMMethods.h"
#include "IDipolarField.h"

class DipolarInteractions : public NMMethods, public IDipolarField {
public:
    DipolarInteractions(std::shared_ptr<SystemDataContainer> data) : NMMethods(data) {}
    ~DipolarInteractions() = default;
public:
    // ####################################            Define Public Variables            ###################################

    // ####################################            Define Public Functions            ###################################
        // Description missing
    double dipolarKernel1D(const int& originSite, const int& influencingSite, const std::string& component) override;
    double dipolarKernel3D(const int& originSite, const int& influencingSite, const double& A, const double& alpha) override;
    void DipolarInteraction1D(std::vector<double> inMxTerms, std::vector<double>& outDipoleX) override;

    std::vector<double> DipolarInteractionClassic(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                     std::vector<double> mzTerms, std::vector<int> sitePositions) override;

    // Description missing
    std::vector<double> DipolarInteractionIntralayer(std::vector<std::vector<double>>& mTerms,
                                                     int& currentSite, const int& currentLayer=0,
                                                     const double& exchangeStiffness=5.3e-17) override;

    // Description missing
    std::vector<double> DipolarInteractionInterlayer(std::vector<std::vector<double>>& mTermsLayer1,
                                                     std::vector<std::vector<double>>& mTermsLayer2, int& currentSite,
                                                     const int& currentLayer, const int& otherLayer) override;

    // Description missing
    std::vector<double> DipolarInteractionInterlayerAll(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2,
                                                        int& currentSite, const int& currentLayer, const int& otherLayer,
                                                        double& exchangeStiffness, double& interlayerExchange) override;

    // Description missing
    std::vector<double> DipolarInteractionInterlayerAdjacent(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2, int& numNeighbours,
                                                        int& currentSite, const int& currentLayer,
                                                        double& exchangeStiffness, double& interlayerExchange) override;

    std::vector<double> DipolarInteractionIntralayerDebug(std::vector<std::vector<double>>& mTerms, int& numNeighbours,
                                              int& currentSite, const int& currentLayer = 0) override;
    // Description missing
    std::vector<double> DipolarInteractionInterlayerDebug(std::vector<std::vector<double>>& mTermsChain1,
                                                     std::vector<std::vector<double>>& mTermsChain2, int& numNeighbours,
                                                     int& currentSite, const int& currentLayer = 0) override;
};


#endif //SPINCHAINS_DIPOLARFIELD_H
