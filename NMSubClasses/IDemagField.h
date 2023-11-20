//
// Created by Cameron Aidan McEleney on 19/11/2023.
//

#ifndef SPINCHAINS_IDEMAGFIELD_H
#define SPINCHAINS_IDEMAGFIELD_H

#include <vector>

class IDemagField {
public:
    ~IDemagField() = default;
public:
    virtual void                   DemagField1DComplex(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage) = 0;

    virtual void                    DemagField1DReal(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage) = 0;

    /*
    virtual void                DemagnetisationFieldFFT (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                                 const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                                 const std::vector<double>& mzTerms) = 0;
                                                 */

    virtual void                DemagnetisationFieldIntense (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                              const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                              const std::vector<double>& mzTerms) = 0;


    virtual void                    DemagFieldsUsingDipoles(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                  std::vector<double> mzTerms, std::vector<int> sitePositions,
                                                  std::vector<double>& outDemagX, std::vector<double>& outDemagY,
                                                  std::vector<double>& outDemagZ) = 0;
};

#endif //SPINCHAINS_IDEMAGFIELD_H
