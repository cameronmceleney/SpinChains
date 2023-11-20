//
// Created by Cameron Aidan McEleney on 20/11/2023.
//

#ifndef SPINCHAINS_IDIPOLARFIELD_H
#define SPINCHAINS_IDIPOLARFIELD_H

class IDipolarField {
public:
    // ####################################            Define Public Variables            ###################################

    // ####################################            Define Public Functions            ###################################
        // Description missing
    virtual double dipolarKernel1D(const int& originSite, const int& influencingSite, const std::string& component) = 0;
    virtual double dipolarKernel3D(const int& originSite, const int& influencingSite, const double& A, const double& alpha) = 0;
    virtual void DipolarInteraction1D(std::vector<double> inMxTerms, std::vector<double>& outDipoleX) = 0;

    virtual std::vector<double> DipolarInteractionClassic(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                     std::vector<double> mzTerms, std::vector<int> sitePositions) = 0;

    // Description missing
    virtual std::vector<double> DipolarInteractionIntralayer(std::vector<std::vector<double>>& mTerms,
                                                     int& currentSite, const int& currentLayer=0,
                                                     const double& exchangeStiffness=5.3e-17) = 0;

    // Description missing
    virtual std::vector<double> DipolarInteractionInterlayer(std::vector<std::vector<double>>& mTermsLayer1,
                                                     std::vector<std::vector<double>>& mTermsLayer2, int& currentSite,
                                                     const int& currentLayer, const int& otherLayer) = 0;

    // Description missing
    virtual std::vector<double> DipolarInteractionInterlayerAll(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2,
                                                        int& currentSite, const int& currentLayer, const int& otherLayer,
                                                        double& exchangeStiffness, double& interlayerExchange) = 0;

    // Description missing
    virtual std::vector<double> DipolarInteractionInterlayerAdjacent(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2, int& numNeighbours,
                                                        int& currentSite, const int& currentLayer,
                                                        double& exchangeStiffness, double& interlayerExchange) = 0;

    virtual std::vector<double> DipolarInteractionIntralayerDebug(std::vector<std::vector<double>>& mTerms, int& numNeighbours,
                                              int& currentSite, const int& currentLayer = 0) = 0;
    // Description missing
    virtual std::vector<double> DipolarInteractionInterlayerDebug(std::vector<std::vector<double>>& mTermsChain1,
                                                     std::vector<std::vector<double>>& mTermsChain2, int& numNeighbours,
                                                     int& currentSite, const int& currentLayer = 0) = 0;

    ~IDipolarField() = default;
};

#endif //SPINCHAINS_IDIPOLARFIELD_H
