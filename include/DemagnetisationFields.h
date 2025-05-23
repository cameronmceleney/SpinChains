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

    enum class DemagMethod {
        DipoleGreenFunctionReal,
        DipoleGreenFunctionComplex,
        DipoleBruteForce,
        Demag
    };
public:
    void initialise(const bool &useDebug);

    template <typename T>
    void
    calculateOneDimension( const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                                  const std::vector<double> &mzTerms, std::vector<T> &demagFieldXOut,
                                                  std::vector<T> &demagFieldYOut, std::vector<T> &demagFieldZOut,
                                                  const DemagnetisationFields::DemagMethod &demagMethod,
                                                  CommonStructures::Parallelisations parallelFlag );


    CommonStructures::Vector3D calculateDemagField(int dimension, int site,
                                                                 const std::vector<double>& mxTerms,
                                                                 const std::vector<double>& myTerms,
                                                                 const std::vector<double>& mzTerms);

    static CommonStructures::Vector3D calculateDipoleFieldBruteForce(int dimension, int site,
                                                             const std::vector<double>& mxTerms,
                                                             const std::vector<double>& myTerms,
                                                             const std::vector<double>& mzTerms);

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    calculateDipoleFieldGreenReal( int dimension, int site,
                                   const std::vector<double>& mxTerms,
                                   const std::vector<double>& myTerms,
                                   const std::vector<double>& mzTerms);

    CommonStructures::Vector3D calculateDipoleFieldGreenComplex( int dimension, int site,
                                                     const std::vector<double>& mxTerms,
                                                     const std::vector<double>& myTerms,
                                                     const std::vector<double>& mzTerms);


private:
    CommonStructures::Vector3D _calculateFieldForSite(const int &site, const std::vector<double>& mxTerms,
                                                                    const std::vector<double>& myTerms,
                                                                    const std::vector<double>& mzTerms,
                                                                    const DemagnetisationFields::DemagMethod& demagMethod);
    static float _calculateDemagFactorsUniformPrism( const double &a, const double &b, const double &c);

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    _calculateDipoleFieldGreen1DReal( const int &site, const std::vector<double> &inMxTerms,
                                      const std::vector<double> &inMyTerms,
                                      const std::vector<double> &inMzTerms );

    CommonStructures::Vector3D _calculateDemagField(int site, const std::vector<double>& mxTerms,
                                               const std::vector<double>& myTerms,
                                               const std::vector<double>& mzTerms);

    static CommonStructures::Vector3D _calculateDipoleFieldBruteForce(const int &site, const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                                          const std::vector<double> &mzTerms);

    CommonStructures::Vector3D _calculateDipoleFieldGreen1DComplex( const int &site, const std::vector<double> &inMxTerms,
                                                               const std::vector<double> &inMyTerms,
                                                               const std::vector<double> &inMzTerms );

    template <typename T>
    void _addDemagField(int site, const CommonStructures::Vector3D& demagField, std::vector<T>& demagFieldXOut,
                        std::vector<T>& demagFieldYOut, std::vector<T>& demagFieldZOut);


    // Template function declaration for memory allocation and checking
    template <typename T>
    T* fftw_alloc_and_check(const char* var_name, int size);

    // Function to populate memory for magnetization components
    void populateMemory(double* mX, double* mY, double* mZ, const std::vector<double>& inMxTerms,
                        const std::vector<double>& inMyTerms, const std::vector<double>& inMzTerms,
                        bool &shouldSkipBounds);
    void populateMemory(fftw_complex* mX, fftw_complex* mY, fftw_complex* mZ, const std::vector<double>& inMxTerms,
                        const std::vector<double>& inMyTerms, const std::vector<double>& inMzTerms,
                        bool& shouldSkipBounds);

    // Function to create and execute FFTW plans
    void createAndExecuteFFTWPlans( fftw_complex *planXIn, fftw_complex *planYIn, fftw_complex *planZIn,
                                    fftw_complex (*planXOut), fftw_complex (*planYOut),
                                    fftw_complex (*planZOut), int trueNumSpins, bool forward );

    // Function to compute Green's function in Fourier space
    void computeGreenFunctionInFourierSpace( double* Gxx, double* Gyy, double* Gzz, int trueNumSpins);
    void computeGreenFunctionInFourierSpace( fftw_complex * Gxx, fftw_complex* Gyy, fftw_complex* Gzz, int trueNumSpins);

    // Function to perform multiplication in Fourier space
    void multiplyInFourierSpace( double* hdX, double* hdY, double* hdZ,
                                 const double* Gxx, const double* Gyy, const double* Gzz,
                                 const fftw_complex *mX, const fftw_complex *mY, const fftw_complex *mZ,
                                 int trueNumSpins);
    void multiplyInFourierSpace( fftw_complex * hdX, fftw_complex* hdY, fftw_complex* hdZ, const fftw_complex* Gxx, const fftw_complex* Gyy, const fftw_complex* Gzz, const fftw_complex* mX, const fftw_complex* mY, const fftw_complex* mZ, int trueNumSpins);

    // Function to compute RMSE between the computed and input magnetization terms
    void computeRMSE(double& rmse_mx, double& rmse_my, double& rmse_mz, const double* mX, const double* mY, const double* mZ, const std::vector<double>& inMxTerms, const std::vector<double>& inMyTerms, const std::vector<double>& inMzTerms, int gotNumSpins);

    // Function to populate the output arrays with demagnetization field components
    void populateOutputArrays(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ, const double* hdX, const double* hdY, const double* hdZ, bool applyDemag, int gotNumSpins, int trueNumSpins);
};


#endif //SPINCHAINS_DEMAGNETISATIONFIELDS_H
