//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/SolversConfiguration.h"

SolversConfiguration::SolversConfiguration(std::shared_ptr<SimulationParameters> sharedSimParams,
                                   std::shared_ptr<SimulationStates> sharedSimStates,
                                   std::shared_ptr<SimulationFlags> sharedSimFlags)

    : SolversSuperClass(std::move(sharedSimParams), std::move(sharedSimStates), std::move(sharedSimFlags)) {
    _mxInit = 0.0;
    _myInit = 0.0;
    _mzInit = 1.0;
    _zeroValue = 0.0;
}

void SolversConfiguration::Configure() {


    // ###################### Core Method Invocations ######################
    // Order is intentional, and must be maintained!
    _testShockwaveInitConditions();

    if (simFlags->hasMultipleLayers) {
        _generateMultilayerAbsorbingRegions(simParams->numSpinsInABC, simParams->gilbertDamping,
                                            simParams->gilbertABCInner, simParams->gilbertABCOuter);
    } else {
        _generateAbsorbingRegions(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->gilbertDamping,
                                  simParams->gilbertABCInner, simParams->gilbertABCOuter);
    }

    _setupDrivingRegion(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->drivingRegionWidth);
    _generateExchangeVector(simParams->numSpinsInABC, simParams->numberOfSpinPairs, simParams->exchangeEnergyMin,
                            simParams->exchangeEnergyMax);

    if (simFlags->hasMultipleLayers) {
        simStates->m0Nest = InitialiseNestedVectors(simParams->numLayers, _mxInit, _myInit, _mzInit);
        simStates->m1Nest = InitialiseNestedVectors(simParams->numLayers, _mxInit, _myInit, _zeroValue);
        simStates->m2Nest = InitialiseNestedVectors(simParams->numLayers, _mxInit, _myInit, _zeroValue);
    } else if (!simFlags->hasMultipleLayers) {
        _setupInitMagneticMoments(_mxInit, _myInit, _mzInit);
    }
}

void SolversConfiguration::_generateAbsorbingRegions(int numSpinsInChain, int numSpinsAbsorbingRegion, double gilbertSpinChain,
                                                     double gilbertAbsorbingRegionInner, double gilbertAbsorbingRegionOuter) {
    // Generate the damping regions that are appended to either end of the spin chain.

    // These damping regions are part of my Periodic Boundary Conditions (PBCs) as Absorbing Boundary Conditions (ABCs).
    LinspaceClass DampingRegionLeft;
    LinspaceClass DampingRegionRight;

    if (numSpinsAbsorbingRegion < 0) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(1);
    }

    std::vector<double> gilbertChain(numSpinsInChain, gilbertSpinChain);

    DampingRegionLeft.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion, true, false);
    DampingRegionRight.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion, true, false);
    std::vector<double> dampingRegionLHS = DampingRegionLeft.generate_array();
    std::vector<double> dampingRegionRHS = DampingRegionRight.generate_array();

    // Combine all damped regions to form vector which describes the entire spinchain.
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionLHS.begin(), dampingRegionLHS.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionRHS.begin(), dampingRegionRHS.end());
    simStates->gilbertVector.push_back(0);
}

void SolversConfiguration::_setupDrivingRegion(int numSpinsInChain, int numSpinsAbsorbingRegion, int drivingRegionWidth) {

    if (simFlags->shouldDriveDiscreteSites) {
        for (int& drivenSite: simStates->discreteDrivenSites)
            drivenSite -= 1;

        // Update so that output file is intuitive
        simParams->drivingRegionWidth = 0;
        simParams->drivingRegionLhs = 0;
        simParams->drivingRegionRhs = 0;
        return;
    }

    if (simFlags->shouldDriveCentre) {
        if (drivingRegionWidth % 2 == 0) {
            simParams->drivingRegionLhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion - (drivingRegionWidth / 2);
            simParams->drivingRegionRhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion + (drivingRegionWidth / 2);
        } else {
            // Use convention to round down then add additional site (odd number) to the RHS of the centre
            simParams->drivingRegionLhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion - (drivingRegionWidth / 2);
            simParams->drivingRegionRhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion + (1 + drivingRegionWidth / 2);
        }
        return;
    }

    if (simFlags->shouldDriveBothSides) {throw std::runtime_error(std::string("shouldDriveBothSides has been selected, but this feature is not yet implemented"));}

    if (simFlags->shouldDriveLHS) {
        // The +1/-1 offset excludes the zeroth spin while retaining the correct driving width
        simParams->drivingRegionLhs = numSpinsAbsorbingRegion + 1;
        simParams->drivingRegionRhs = simParams->drivingRegionLhs + drivingRegionWidth - 1;
        return;
    }

    if (simFlags->shouldDriveRHS) {
        // The +1 is to correct the offset of adding a zeroth spin
        simParams->drivingRegionRhs = simParams->systemTotalSpins - numSpinsAbsorbingRegion - 1;
        simParams->drivingRegionLhs = simParams->drivingRegionRhs - drivingRegionWidth + 1;
        return;
    }

    if (simFlags->hasCustomDrivePosition) {
        // The +1/-1 offset excludes the zeroth spin while retaining the correct driving width
        simParams->drivingRegionLhs = simParams->numSpinsInABC + 1;
        simParams->drivingRegionRhs = simParams->drivingRegionLhs + drivingRegionWidth - 1;
        return;
    }

}
void SolversConfiguration::_generateExchangeVector(int numSpinsAbsorbingRegion, int numSpinPairs, double exchangeMin, double exchangeMax) {
    /*
     * Create the arrays which house the exchange integral values. There are options to have a non-uniform exchange coded in, as well as the option to
     * induce a 'kick' into the system by initialising certain spins to have differing parameters to their neighbours.
     */
    LinspaceClass SpinChainExchange;
    LinspaceClass SpinChainExchangeLHS;
    LinspaceClass SpinChainExchangeRHS;

    int customDampedSite = -1; // Lets the damping regions be a single exchange region with the main chain

    if (numSpinsAbsorbingRegion > 0) {
        if (simParams->numSpinsInABC == customDampedSite) {
            numSpinPairs += (2 * simParams->numSpinsInABC);
            SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, true);
            simStates->exchangeVec = SpinChainExchange.generate_array();
        } else {
            simStates->exchangeVec.push_back(0.0);

            SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, false);
            std::vector<double> tempExchangeChain = SpinChainExchange.generate_array();

            SpinChainExchangeLHS.set_values(exchangeMin, exchangeMax, numSpinsAbsorbingRegion, true, false);
            std::vector<double> tempExchangeLHS = SpinChainExchangeLHS.generate_array();

            SpinChainExchangeRHS.set_values(exchangeMax, exchangeMin, numSpinsAbsorbingRegion, true, false);
            std::vector<double> tempExchangeRHS = SpinChainExchangeRHS.generate_array();

            simStates->exchangeVec.insert(simStates->exchangeVec.end(), tempExchangeLHS.begin(), tempExchangeLHS.end());
            simStates->exchangeVec.insert(simStates->exchangeVec.end(), tempExchangeChain.begin(), tempExchangeChain.end());
            simStates->exchangeVec.insert(simStates->exchangeVec.end(), tempExchangeRHS.begin(), tempExchangeRHS.end());
            simStates->exchangeVec.push_back(0.0);
        }
    } else {
        // The linearly spaced vector is saved as the class member 'exchangeVec' simply to increase code readability
        SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, true);
        simStates->exchangeVec = SpinChainExchange.generate_array();
    }
}
void SolversConfiguration::_setupInitMagneticMoments(double mxInit, double myInit, double mzInit) {

    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mxInitCond(simParams->systemTotalSpins, mxInit), myInitCond(simParams->systemTotalSpins, myInit), mzInitCond(simParams->systemTotalSpins, mzInit);
    // mxInitCond[0] = mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < systemTotalSpins; i++) {
        InitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    if (!simFlags->isFerromagnetic) {
        for (int i = 0; i < simParams->systemTotalSpins; i++) {
            if (i % 2 == 1)
                mzInitCond[i] *= -1.0;
        }
    }

    // Appends initial conditions to the vectors
    simStates->mx0.insert(simStates->mx0.end(), mxInitCond.begin(), mxInitCond.end());
    simStates->my0.insert(simStates->my0.end(), myInitCond.begin(), myInitCond.end());
    simStates->mz0.insert(simStates->mz0.end(), mzInitCond.begin(), mzInitCond.end());

    // This zero is the (N+1)th spin on the RHS of the chain
    simStates->mx0.push_back(0);
    simStates->my0.push_back(0);
    simStates->mz0.push_back(0);
}
void SolversConfiguration::_testShockwaveInitConditions() {

    if (simFlags->hasShockwave) {
        simParams->shockwaveStepsize = (simParams->shockwaveMax - simParams->shockwaveInitialStrength) / simParams->shockwaveGradientTime;
    } else {
        // Ensures, on the output file, all parameter read as zero; reduces confusion when no shockwave is applied.
        simParams->iterStartShock = 0.0;
        simParams->shockwaveScaling = 0.0;
        simParams->shockwaveGradientTime = 0.0;
        simParams->shockwaveInitialStrength = 0.0;
        simParams->shockwaveMax = 0.0;
        simParams->shockwaveStepsize = 0.0;
    }
}

void SolversConfiguration::_generateMultilayerAbsorbingRegions(int numSpinsAbsorbingRegion, double gilbertSpinChain,
                                                          double gilbertAbsorbingRegionInner, double gilbertAbsorbingRegionOuter) {

    LinspaceClass AbsorbingRegionLHS;
    LinspaceClass AbsorbingRegionRHS;

    if (numSpinsAbsorbingRegion < 0) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(1);
    }

    for (int& layer: simStates->layerSpinsInChain) {
        std::vector<double> gilbertChain(layer, gilbertSpinChain);

        AbsorbingRegionLHS.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion, true, false);
        AbsorbingRegionRHS.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion, true, false);

        // Change of name to be specific that this is only the damping region
        std::vector<double> dampingRegionLHS = AbsorbingRegionLHS.generate_array();
        std::vector<double> dampingRegionRHS = AbsorbingRegionRHS.generate_array();

        // Combine all damped regions to form vector which describes the entire spinchain.
        simStates->gilbertVectorMulti[layer].insert(simStates->gilbertVectorMulti[layer].end(), dampingRegionLHS.begin(), dampingRegionLHS.end());
        simStates->gilbertVectorMulti[layer].insert(simStates->gilbertVectorMulti[layer].end(), gilbertChain.begin(), gilbertChain.end());
        simStates->gilbertVectorMulti[layer].insert(simStates->gilbertVectorMulti[layer].end(), dampingRegionRHS.begin(), dampingRegionRHS.end());
        simStates->gilbertVectorMulti[layer].push_back(0);

        //PrintVector(gilbertVectorMulti[i], false);
    }
}
void SolversConfiguration::_setupInitMultilayerMagneticMoments(std::vector<std::vector<std::vector<double>>>& nestedNestedVector,
                                                          int layer, double mxInit, double myInit, double mzInit) {

    // mxInitCond[0] = mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < systemTotalSpins; i++) {
        mxInitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    for (int i = 0; i < simStates->layerTotalSpins[layer]; i++) {
        // i is simply an index here
        nestedNestedVector[layer].push_back({mxInit, myInit, mzInit});
    }

    // This zero is the (N+1)th spin on the RHS of the chain
    nestedNestedVector[layer].push_back({0.0, 0.0, 0.0});

}
std::vector<std::vector<std::vector<double>>> SolversConfiguration::initializeNestedNestedVector(int numSites,
                                                                                                 bool includeEnd) {
    /* Legacy code not used currently. Example implementation is below
     * std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested3;
     * mValsNested3["nestedNestedVector3"] = initializeNestedNestedVector(1, true);
     * std::vector<std::vector<std::vector<double>>> m2Nest = mValsNested3["nestedNestedVector3"];
     * _setupInitMultilayerMagneticMoments(m2Nest, 1, 0, 0 , 0);
     */
    std::vector<std::vector<std::vector<double>>> innerNestedVector;
    for (int j = 0; j < numSites; j++) {
        std::vector<std::vector<double>> innerVector;
        std::vector<double> innermostVector = {0.0, 0.0, 0.0};
        innerVector.push_back(innermostVector);
        innerNestedVector.push_back(innerVector);
    }
    return innerNestedVector;
}

std::vector<std::vector<std::vector<double>>> SolversConfiguration::InitialiseNestedVectors(int& totalLayer, double& mxInit, double& myInit, double& mzInit) {

    // Initialise mapping
    std::map<std::string, std::vector<std::vector<std::vector<double>>>> mTermsMapping;

    // This is likely a very slow way to initialise (push_back is slow), but this works for now. Fix if it is a bottleneck later
    std::vector<std::vector<std::vector<double>>> innerNestedVector;
    for (int i = 0; i < totalLayer; i++) {
        std::vector<std::vector<double>> innerVector;
        std::vector<double> innermostVector = {0.0, 0.0, 0.0};
        innerVector.push_back(innermostVector);
        innerNestedVector.push_back(innerVector);
    }
    // Assign name to nested-nested vector
    mTermsMapping["nestedVector"] = innerNestedVector;

    // Assign key of map to multi-dim vector
    std::vector<std::vector<std::vector<double>>> mTermsNested = mTermsMapping["nestedVector"];

    // Invoke method to set initial magnetic moments. To call: mValsNest[layer][site][component]
    for (int layer = 0; layer < totalLayer; layer++)
        _setupInitMultilayerMagneticMoments(mTermsNested, layer, mxInit, myInit , mzInit);

    return mTermsNested;
}