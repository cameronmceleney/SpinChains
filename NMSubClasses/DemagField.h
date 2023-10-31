//
// Created by Cameron McEleney on 31/10/2023.
//
#pragma once

#ifndef SPINCHAINS_DEMAGFIELD_H
#define SPINCHAINS_DEMAGFIELD_H

#include "NMMethods.h"

class DemagnetisationFields:  public NMMethods {
private:
    // ####################################            Define Private Variables            ###################################

    // ####################################            Define Private Methods            ###################################

protected:
    // ####################################            Define Protected Variables            ###################################

    // ####################################            Define Protected Methods            ###################################


public:
    // ####################################            Define Public Variables            ###################################

    // ####################################            Define Public Methods            ###################################
    void                    DemagField1DComplex(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage);

    void                    DemagField1DReal(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage);

    void                DemagnetisationFieldFFT (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                                 const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                                 const std::vector<double>& mzTerms);

    void                DemagnetisationFieldIntense (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                              const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                              const std::vector<double>& mzTerms);

    void                    DemagFieldsUsingDipoles(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                  std::vector<double> mzTerms, std::vector<int> sitePositions,
                                                  std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ);

};


#endif //SPINCHAINS_DEMAGFIELD_H
