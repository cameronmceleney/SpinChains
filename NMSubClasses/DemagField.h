//
// Created by Cameron McEleney on 31/10/2023.
//
#pragma once

#ifndef SPINCHAINS_DEMAGFIELD_H
#define SPINCHAINS_DEMAGFIELD_H

#include "NMMethods.h"
#include "IDemagField.h"

class DemagnetisationFields : public NMMethods, public IDemagField {
public:
    DemagnetisationFields(std::shared_ptr<SystemDataContainer> data) : NMMethods(data) {}
public:
    void                    DemagField1DComplex(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage) override;

    void                    DemagField1DReal(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage) override;

    /*
    void                DemagnetisationFieldFFT (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                                 const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                                 const std::vector<double>& mzTerms) override;
                                                 */

    void                DemagnetisationFieldIntense (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                              const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                              const std::vector<double>& mzTerms) override;


    void                    DemagFieldsUsingDipoles(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                  std::vector<double> mzTerms, std::vector<int> sitePositions,
                                                  std::vector<double>& outDemagX, std::vector<double>& outDemagY,
                                                  std::vector<double>& outDemagZ) override;


};


#endif //SPINCHAINS_DEMAGFIELD_H
