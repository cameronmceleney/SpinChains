//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSCONFIGURATION_H
#define SPINCHAINS_SOLVERSCONFIGURATION_H

// C++ Standard Library
#include <iostream>
#include <map>
#include <string>

// C++ User Library (Parent)
#include "SolversSuperClass.h"

// C++ User Library (General)
#include "../libs/linspace.h"

class SolversConfiguration:  public SolversSuperClass {
private:
    // ####################################            Define Private Variables            ###################################

    double _mxInit;  // Initial magnetisation in the x-direction
    double _myInit;  // Initial magnetisation in the y-direction
    double _mzInit;  // Initial magnetisation in the z-direction
    double _zeroValue; // Saves retyping 0.0

    // ####################################            Define Private Functions            ###################################
    // #################################     all systems    ###############################

    /**
     * Set up driving regions for the system. The LHS option is solely for drives from the left of the system. The RHS options contains the
     * drive from the right, as well as an option to drive from the centre.
     **/
    void                _setupDrivingRegion(int numSpinsInChain, int numSpinsAbsorbingRegion, int drivingRegionWidth);

    // Description missing
    void                _generateExchangeVector(int numSpinsAbsorbingRegion, int numSpinPairs, double exchangeMin, double exchangeMax);

    // ##############################      single systems    #############################
    // Description missing
    void                _generateAbsorbingRegions(int numSpinsInChain, int numSpinsAbsorbingRegion, double gilbertSpinChain,
                                                  double gilbertAbsorbingRegionInner, double gilbertAbsorbingRegionOuter);

    // Description missing
    void                _setupInitMagneticMoments(double mxInit, double myInit, double mzInit);

    // ###########################     multilayered systems    ##########################

    // Generate the damping regions (part of the Absorbing Boundary Conditions (ABCs)) that are appended to both ends of the spin chain.
    void                _generateMultilayerAbsorbingRegions(int numSpinsAbsorbingRegion, double gilbertSpinChain,
                                                            double gilbertAbsorbingRegionInner, double gilbertAbsorbingRegionOuter);

    // Description missing
    std::vector<std::vector<std::vector<double>>> initializeNestedNestedVector(int numSites, bool includeEnd);

    // Description missing
    std::vector<std::vector<std::vector<double>>> InitialiseNestedVectors(int& totalLayer, double& mxInit,
                                                                          double& myInit, double& mzInit);


protected:
    void                _testShockwaveInitConditions();

public:
    SolversConfiguration(std::shared_ptr<SimulationParameters> paramsData,
                                   std::shared_ptr<SimulationStates> sharedSimStates,
                                   std::shared_ptr<SimulationFlags> sharedSimFlags);

    ~SolversConfiguration() = default;

public:
    // Assignment of all values required for the simulation
    void                Configure();

    void                _setupInitMultilayerMagneticMoments(std::vector<std::vector<std::vector<double>>>& nestedNestedVector,
                                                            int layer, double mxInit, double myInit, double mzInit);

    void                performInitialisation() override { Configure(); };
};


#endif //SPINCHAINS_SOLVERSCONFIGURATION_H
