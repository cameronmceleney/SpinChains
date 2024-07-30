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
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

// C++ User Libraries (General)
#include "../libs/CommonDefinitions.h"
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

    void                    DemagField1DRealGreen(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                                  std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                                  int iteration, std::string rkStageName);

    /*
    void                    DemagnetisationFieldFFT (std::vector<double>& H_dx, std::vector<double>& H_dy, std::vector<double>& H_dz,
                                                 const std::vector<double>&mxTerms, const std::vector<double>& myTerms,
                                                 const std::vector<double>& mzTerms);
                                                 */

    void                    DemagnetisationFieldIntense( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                      const std::vector<double> &mzTerms, std::vector<double> &H_dx,
                                      std::vector<double> &H_dy, std::vector<double> &H_dz );


    void                    DemagFieldsUsingDipoles(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                    std::vector<double> mzTerms, std::vector<int> sitePositions,
                                                    std::vector<double>& outDemagX, std::vector<double>& outDemagY,
                                                    std::vector<double>& outDemagZ);

    void                    calculateOneDimension( const std::vector<double> &mxTermsIn,
                                                   const std::vector<double> &myTermsIn,
                                                   const std::vector<double> &mzTermsIn,
                                                   std::vector<std::atomic<double>> &demagFieldXOut,
                                                   std::vector<std::atomic<double>> &demagFieldYOut,
                                                   std::vector<std::atomic<double>> &demagFieldZOut,
                                                   const bool &shouldUseTBB );

void                        initialise();

private:
    float _calculateDemagFactorsUniformPrism( const double &a, const double &b, const double &c);
};


#endif //SPINCHAINS_DEMAGNETISATIONFIELDS_H
