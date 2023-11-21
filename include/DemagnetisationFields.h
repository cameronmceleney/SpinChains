//
// Created by Cameron McEleney on 31/10/2023.
//
#pragma once

#ifndef SPINCHAINS_DEMAGNETISATIONFIELDS_H
#define SPINCHAINS_DEMAGNETISATIONFIELDS_H

// C++ Standard Library
#include <cstring>
#include <iostream>
#include <vector>


// C++ Third Party Library
extern "C"{
    #include <fftw3.h>
}

// C++ User Libraries (General)
#include "CommonLibs.h"
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"

class DemagnetisationFields {
private:
    SimulationParameters* _simParams; // Non-owning pointer to SimulationParameters
    SimulationStates* _simStates;
    SimulationFlags* _simFlags;
public:
    explicit DemagnetisationFields(SimulationParameters* sharedSimParams,
                                   SimulationStates* sharedSimStates,
                                   SimulationFlags* sharedSimFlags);
    ~DemagnetisationFields() = default;
public:
    void                    DemagField1DComplex(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage);

    void                    DemagField1DReal(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                     std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                     int iteration, std::string rkStage);

    /*
    void                DemagnetisationFieldFFT (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                                 const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                                 const std::vector<double>& mzTerms);
                                                 */

    void                DemagnetisationFieldIntense (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                              const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                              const std::vector<double>& mzTerms);


    void                    DemagFieldsUsingDipoles(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                  std::vector<double> mzTerms, std::vector<int> sitePositions,
                                                  std::vector<double>& outDemagX, std::vector<double>& outDemagY,
                                                  std::vector<double>& outDemagZ);


};


#endif //SPINCHAINS_DEMAGNETISATIONFIELDS_H
