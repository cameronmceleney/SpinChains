//
// Created by Cameron McEleney on 31/10/2023.
//

#include "NMConfiguration.h"

void NMConfiguration::Configure() {
    _mxInit = 0.0;
    _myInit = 0.0;
    _mzInit = 1.0;
    double zeroValue = 0.0;

    // ###################### Core Method Invocations ######################
    // Order is intentional, and must be maintained!
    _testShockwaveInitConditions();
    if (systemData->useMultilayer) {
        _generateMultilayerAbsorbingRegions(systemData->numSpinsDamped, systemData->gilbertDamping,
                                            systemData->gilbertABCInner, systemData->gilbertABCOuter);
    } else {
        _generateAbsorbingRegions(systemData->numSpinsInChain, systemData->numSpinsDamped, systemData->gilbertDamping,
                                  systemData->gilbertABCInner, systemData->gilbertABCOuter);}
    _setupDrivingRegion(systemData->numSpinsInChain, systemData->numSpinsDamped, systemData->drivingRegionWidth);
    _generateExchangeVector(systemData->numSpinsDamped, systemData->numberOfSpinPairs, systemData->exchangeEnergyMin,
                            systemData->exchangeEnergyMax);

    if (systemData->useMultilayer) {
         systemData->m0Nest = InitialiseNestedVectors(systemData->totalLayers, _mxInit, _myInit, _mzInit);
         systemData->m1Nest = InitialiseNestedVectors(systemData->totalLayers, _mxInit, _myInit, zeroValue);
         systemData->m2Nest = InitialiseNestedVectors(systemData->totalLayers, _mxInit, _myInit, zeroValue);
    } else if (!systemData->useMultilayer) {
        _setupInitMagneticMoments(_mxInit, _myInit, _mzInit);
    }

    std::cout << "ambientTemp in config: " << systemData->ambientTemperature << std::endl;
    std::cout<<"Completed initialisation of config" << std::endl;
}

void NMConfiguration::_generateAbsorbingRegions(int numSpinsInChain, int numSpinsAbsorbingRegion, double gilbertSpinChain,
                                                double gilbertAbsorbingRegionInner, double gilbertAbsorbingRegionOuter) {
    // Generate the damping regions that are appended to either end of the spin chain.

    // These damping regions are part of my Periodic Boundary Conditions (PBCs) as Absorbing Boundary Conditions (ABCs).
    LinspaceClass DampingRegionLeft;
    LinspaceClass DampingRegionRight;

    if (numSpinsAbsorbingRegion < 0) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(0);
    }

    std::vector<double> gilbertChain(numSpinsInChain, gilbertSpinChain);

    DampingRegionLeft.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion, true, false);
    DampingRegionRight.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion, true, false);
    std::vector<double> dampingRegionLHS = DampingRegionLeft.generate_array();
    std::vector<double> dampingRegionRHS = DampingRegionRight.generate_array();

    // Combine all damped regions to form vector which describes the entire spinchain.
    systemData->gilbertVector.insert(systemData->gilbertVector.end(), dampingRegionLHS.begin(), dampingRegionLHS.end());
    systemData->gilbertVector.insert(systemData->gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());
    systemData->gilbertVector.insert(systemData->gilbertVector.end(), dampingRegionRHS.begin(), dampingRegionRHS.end());
    systemData->gilbertVector.push_back(0);
}

void NMConfiguration::_setupDrivingRegion(int numSpinsInChain, int numSpinsAbsorbingRegion, int drivingRegionWidth) {

    if (systemData->centralDrive) {
        systemData->drivingRegionLhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion - (drivingRegionWidth / 2);
        systemData->drivingRegionRhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion + (drivingRegionWidth / 2);
        return;
    }

    if (systemData->dualDrive) {throw std::runtime_error(std::string("dualDrive has been selected, but this feature is not yet implemented"));}

    if (systemData->lhsDrive) {
        // The +1/-1 offset excludes the zeroth spin while retaining the correct driving width
        systemData->drivingRegionLhs = numSpinsAbsorbingRegion + 1;
        systemData->drivingRegionRhs = systemData->drivingRegionLhs + drivingRegionWidth - 1;
        return;
    }

    if (!systemData->rhsDrive) {
        // The +1 is to correct the offset of adding a zeroth spin
        systemData->drivingRegionRhs = systemData->systemTotalSpins - numSpinsAbsorbingRegion - 1;
        systemData->drivingRegionLhs = systemData->drivingRegionRhs - drivingRegionWidth + 1;
        return;
    }

}
void NMConfiguration::_generateExchangeVector(int numSpinsAbsorbingRegion, int numSpinPairs, double exchangeMin, double exchangeMax) {
    /*
     * Create the arrays which house the exchange integral values. There are options to have a non-uniform exchange coded in, as well as the option to
     * induce a 'kick' into the system by initialising certain spins to have differing parameters to their neighbours.
     */
    LinspaceClass SpinChainExchange;

    if (numSpinsAbsorbingRegion > 0) {
        SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, false);
        systemData->exchangeVec = SpinChainExchange.generate_array();

        std::vector<double> dampingRegionLeftExchange(numSpinsAbsorbingRegion, exchangeMin), dampingRegionRightExchange(numSpinsAbsorbingRegion, exchangeMax);
        dampingRegionLeftExchange.insert(dampingRegionLeftExchange.begin(), 0);
        dampingRegionRightExchange.push_back(0);

        systemData->exchangeVec.insert(systemData->exchangeVec.begin(), dampingRegionLeftExchange.begin(), dampingRegionLeftExchange.end());
        systemData->exchangeVec.insert(systemData->exchangeVec.end(), dampingRegionRightExchange.begin(), dampingRegionRightExchange.end());
    } else {
        // The linearly spaced vector is saved as the class member 'exchangeVec' simply to increase code readability
        SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, true);
        systemData->exchangeVec = SpinChainExchange.generate_array();
    }
}
void NMConfiguration::_setupInitMagneticMoments(double mxInit, double myInit, double mzInit) {

    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mxInitCond(systemData->systemTotalSpins, mxInit), myInitCond(systemData->systemTotalSpins, myInit), mzInitCond(systemData->systemTotalSpins, mzInit);
    // mxInitCond[0] = mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < systemTotalSpins; i++) {
        InitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    if (!systemData->isFm) {
        for (int i = 0; i < systemData->systemTotalSpins; i++) {
            if (i % 2 == 1)
                mzInitCond[i] *= -1.0;
        }
    }

    // Appends initial conditions to the vectors
    systemData->mx0.insert(systemData->mx0.end(), mxInitCond.begin(), mxInitCond.end());
    systemData->my0.insert(systemData->my0.end(), myInitCond.begin(), myInitCond.end());
    systemData->mz0.insert(systemData->mz0.end(), mzInitCond.begin(), mzInitCond.end());

    // This zero is the (N+1)th spin on the RHS of the chain
    systemData->mx0.push_back(0);
    systemData->my0.push_back(0);
    systemData->mz0.push_back(0);
}
void NMConfiguration::_testShockwaveInitConditions() {

    if (systemData->hasShockwave) {
        systemData->shockwaveStepsize = (systemData->shockwaveMax - systemData->shockwaveInitialStrength) / systemData->shockwaveGradientTime;
    } else {
        // Ensures, on the output file, all parameter read as zero; reduces confusion when no shockwave is applied.
        systemData->iterStartShock = 0.0;
        systemData->shockwaveScaling = 0.0;
        systemData->shockwaveGradientTime = 0.0;
        systemData->shockwaveInitialStrength = 0.0;
        systemData->shockwaveMax = 0.0;
        systemData->shockwaveStepsize = 0.0;
    }
}

void NMConfiguration::_generateMultilayerAbsorbingRegions(int numSpinsAbsorbingRegion, double gilbertSpinChain,
                                                          double gilbertAbsorbingRegionInner, double gilbertAbsorbingRegionOuter) {


    LinspaceClass AbsorbingRegionLHS;
    LinspaceClass AbsorbingRegionRHS;

    if (numSpinsAbsorbingRegion < 0) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(0);
    }

    for (int& layer: systemData->layerSpinsInChain) {
        std::vector<double> gilbertChain(layer, gilbertSpinChain);

        AbsorbingRegionLHS.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion, true, false);
        AbsorbingRegionRHS.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion, true, false);

        // Change of name to be specific that this is only the damping region
        std::vector<double> dampingRegionLHS = AbsorbingRegionLHS.generate_array();
        std::vector<double> dampingRegionRHS = AbsorbingRegionRHS.generate_array();

        // Combine all damped regions to form vector which describes the entire spinchain.
        systemData->gilbertVectorMulti[layer].insert(systemData->gilbertVectorMulti[layer].end(), dampingRegionLHS.begin(), dampingRegionLHS.end());
        systemData->gilbertVectorMulti[layer].insert(systemData->gilbertVectorMulti[layer].end(), gilbertChain.begin(), gilbertChain.end());
        systemData->gilbertVectorMulti[layer].insert(systemData->gilbertVectorMulti[layer].end(), dampingRegionRHS.begin(), dampingRegionRHS.end());
        systemData->gilbertVectorMulti[layer].push_back(0);

        //PrintVector(gilbertVectorMulti[i], false);
    }
}
void NMConfiguration::_setupInitMultilayerMagneticMoments(std::vector<std::vector<std::vector<double>>>& nestedNestedVector,
                                                          int layer, double mxInit, double myInit, double mzInit) {

    // mxInitCond[0] = mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < systemTotalSpins; i++) {
        mxInitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    for (int i = 0; i < systemData->layerTotalSpins[layer]; i++) {
        // i is simply an index here
        nestedNestedVector[layer].push_back({mxInit, myInit, mzInit});
    }

    // This zero is the (N+1)th spin on the RHS of the chain
    nestedNestedVector[layer].push_back({0.0, 0.0, 0.0});

}

std::vector<std::vector<std::vector<double>>> NMConfiguration::InitialiseNestedVectors(int& totalLayer, double& mxInit, double& myInit, double& mzInit) {

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
std::vector<std::vector<std::vector<double>>> NMConfiguration::initializeNestedNestedVector(int numLayers, bool includeEnd) {
    /* Legacy code, not used in current implementation. Example implementation is below

    std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested3;
    mValsNested3["nestedNestedVector3"] = initializeNestedNestedVector(1, true);
    std::vector<std::vector<std::vector<double>>> m2Nest = mValsNested3["nestedNestedVector3"];
    _setupInitMultilayerMagneticMoments(m2Nest, 1, 0, 0 , 0);
    */
    std::vector<std::vector<std::vector<double>>> innerNestedVector;
    for (int j = 0; j < numLayers; j++) {
        std::vector<std::vector<double>> innerVector;
        std::vector<double> innermostVector = {0.0, 0.0, 0.0};
        innerVector.push_back(innermostVector);
        innerNestedVector.push_back(innerVector);
    }
    return innerNestedVector;
}