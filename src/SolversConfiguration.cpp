//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/SolversConfiguration.h"

SolversConfiguration::SolversConfiguration( std::shared_ptr<SimulationManager> sharedSimManager,
                                            std::shared_ptr<SimulationParameters> sharedSimParams,
                                            std::shared_ptr<SimulationStates> sharedSimStates,
                                            std::shared_ptr<SimulationFlags> sharedSimFlags )

        : SolversSuperClass(std::move(sharedSimManager), std::move(sharedSimParams), std::move(sharedSimStates),
                            std::move(sharedSimFlags)) {
    _mxInit = 0.0;
    _myInit = 0.0;
    _mzInit = 1.0;
    _zeroValue = 0.0;
}

void SolversConfiguration::configure() {

    // ###################### Core Method Invocations ######################
    // Order is intentional, and must be maintained!

    if ( simManager->hasFirstRunOccurred ) {
        // Idiot-proofing from myself
        _reconfigureSystem();
    }

    _testShockwaveInitConditions();

    _setupDrivingRegion(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->drivingRegionWidth);
    _generateExchangeVector(simParams->numSpinsInABC, simParams->numberOfSpinPairs, simParams->exchangeEnergyMin,
                            simParams->exchangeEnergyMax);

    if ( simFlags->hasMultipleStackedLayers ) {
        _generateMultilayerAbsorbingRegions(simParams->numSpinsInABC, simParams->gilbertDamping,
                                            simParams->gilbertABCInner, simParams->gilbertABCOuter);

        simStates->m0Nest = InitialiseNestedVectors(simParams->numLayers, _mxInit, _myInit, _mzInit);
        simStates->m1Nest = InitialiseNestedVectors(simParams->numLayers, _mxInit, _myInit, _zeroValue);
        simStates->m2Nest = InitialiseNestedVectors(simParams->numLayers, _mxInit, _myInit, _zeroValue);
    } else {
        if (simFlags->hasGradientRegionForOscillatingZeeman) {
            if (simFlags->useGenerateABCUpdated)
                _generateAbsorbingRegionsUpdated(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->gilbertDamping,
                                                 simParams->gilbertABCInner, simParams->gilbertABCOuter);
            else
                _generateAbsorbingRegions(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->numSpinsDRPeak,
                                          simParams->numSpinsDRGradient, simParams->gilbertDamping,
                                          simParams->dampingGradientPeak, simParams->gilbertABCInner, simParams->gilbertABCOuter);
        }
        else {
            //_generateAbsorbingRegions(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->gilbertDamping,
             //                         simParams->gilbertABCInner, simParams->gilbertABCOuter);
            _generateAbsorbingRegionsUpdated(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->gilbertDamping,
                                                 simParams->gilbertABCInner, simParams->gilbertABCOuter);
        }

        _setupInitMagneticMoments(_mxInit, _myInit, _mzInit);
        _generateMaps();
    }

    //PrintVector(simStates->dampingGradientMap, true);

    simManager->hasFirstRunOccurred = true;
}

void SolversConfiguration::PrintVector( std::vector<double> &vectorToPrint, bool shouldExitAfterPrint ) {

    std::cout << "\n";

    int count = 0;
    for ( double i: vectorToPrint ) {
        if ( ++count % 10 == 0 )
            std::cout << std::setw(8) << i << std::endl;
        else
            std::cout << std::setw(8) << i << ", ";
    }
    std::cout << "\n\n";

    if ( shouldExitAfterPrint )
        exit(0);
}
void SolversConfiguration::PrintVector( std::vector<int> &vectorToPrint, bool shouldExitAfterPrint ) {

    std::cout << "\n";

    int count = 0;
    for ( double i: vectorToPrint ) {
        if ( ++count % 10 == 0 )
            std::cout << std::setw(8) << i << std::endl;
        else
            std::cout << std::setw(8) << i << ", ";
    }
    std::cout << "\n\n";

    if ( shouldExitAfterPrint )
        exit(0);
}
void SolversConfiguration::PrintVector( std::map<int, double> &mapToPrint, bool shouldExitAfterPrint ) {

    std::cout << "\n";

    int count = 0;
    for ( auto const& [key, val] : mapToPrint ) {
        if ( ++count % 10 == 0 )
            std::cout << std::setw(4) << key << ": " << std::setw(8) << val << std::endl;
        else
            std::cout << std::setw(4) << key << ": " << std::setw(4) << val << ", ";
    }
    std::cout << "\n\n";

    if ( shouldExitAfterPrint )
        exit(0);
}

void SolversConfiguration::_reconfigureSystem() {
    // ###################### Core Method Invocations ######################

    if ( !simManager->hasFirstRunOccurred or simFlags->resetSimState ) {
        // Idiot-proofing from myself
        configure();
    }

    resetInitMagneticMoments(_mxInit, _myInit, _mzInit);
    _setupDrivingRegion(simParams->numSpinsInChain, simParams->numSpinsInABC, simParams->drivingRegionWidth);
}

void SolversConfiguration::_generateMaps() {
    if (simFlags->hasGradientRegionForOscillatingZeeman) {

        // TODO. This sort of implementation is temporary until I convert over to using STRUCTS
        double bias_min, bias_max;
        if (simFlags->isOscillatingZeemanLinearAcrossMap) {
            bias_min = 1.0;
            bias_max = 1.0;
        } else {
            bias_min = 0.0;
            bias_max = 1.0;
        }
        _generateMapOfScaling("Driving Region", simStates->dRGradientMap, simParams->drivingRegionLhs,
                              simParams->drivingRegionRhs, simParams->drivingRegionWidth,
                              simParams->numSpinsDRGradient, simParams->numSpinsDRPeak, bias_min, bias_max);
    }

    if ((simFlags->shouldDmiGradientMirrorOscillatingZeeman || simFlags->shouldDampingGradientMirrorOscillatingZeeman)
    && simStates->dRGradientMap.empty()) {
        throw std::runtime_error("The gradient map is empty, but other simulation flags indicate that it should not be.");
    }

    if (simFlags->shouldDmiGradientMirrorOscillatingZeeman) {
        simParams->dmiRegionLhs = simParams->drivingRegionLhs;
        simParams->dmiRegionRhs = simParams->drivingRegionRhs;
        simParams->numSpinsDmiPeak = simParams->numSpinsDRPeak;
        simParams->numSpinsDmiGradient = simParams->numSpinsDRGradient;
        simParams->numSpinsDmiWidth = simParams->drivingRegionWidth;
        simStates->dmiGradientMap = simStates->dRGradientMap;
    }
    else if ( simFlags->hasGradientRegionForDmi ) {
        if (simParams->dmiRegionLhs < 0 && simParams->dmiRegionRhs < 0) {
            simParams->dmiRegionLhs = simParams->drivingRegionLhs - simParams->dmiRegionOffset;
            simParams->dmiRegionRhs = simParams->drivingRegionRhs + simParams->dmiRegionOffset;
        }
        if ( simParams->numSpinsDmiWidth < 0 ) {
            simParams->numSpinsDmiWidth = simParams->dmiRegionRhs - simParams->dmiRegionLhs + 1;
        }
        simParams->numSpinsDmiGradient = (simParams->numSpinsDmiWidth - simParams->numSpinsDmiPeak) / 2;

        double dmi_min, dmi_max;
        if (simFlags->isDmiLinearAcrossMap) {
            dmi_min = 1.0;
            dmi_max = 1.0;
        } else {
            dmi_min = 0.0;
            dmi_max = 1.0;
        }

        _generateMapOfScaling("DMI Region", simStates->dmiGradientMap, simParams->dmiRegionLhs,
                              simParams->dmiRegionRhs, simParams->numSpinsDmiWidth,
                              simParams->numSpinsDmiGradient, simParams->numSpinsDmiPeak, dmi_min, dmi_max);
    }

    if (simFlags->shouldDampingGradientMirrorOscillatingZeeman) {
        simParams->dampingRegionLhs = simParams->drivingRegionLhs;
        simParams->dampingRegionRhs = simParams->drivingRegionRhs;
        simParams->numSpinsDampingPeak = simParams->numSpinsDRPeak;
        simParams->numSpinsDampingGradient = simParams->numSpinsDRGradient;
        simParams->numSpinsDampingWidth = simParams->drivingRegionWidth;
        simStates->dampingGradientMap = simStates->dRGradientMap;

        if ( !simFlags->isOscillatingZeemanLinearAcrossMap ) {
            throw("SolverConfiguration::_generateMaps: Oscillating Zeeman is not linear across map, but damping gradient "
                  "is mirrored, so this creates a logic error.");
        }

        for (auto const& [key, val] : simStates->dRGradientMap) {
            simStates->dampingGradientMap[key] *= simParams->gilbertDamping;
        }
    }
    else if ( simFlags->hasGradientRegionForDamping ) {
        if (simParams->dampingRegionLhs < 0 && simParams->dampingRegionRhs < 0) {
            simParams->dampingRegionLhs = simParams->drivingRegionLhs - simParams->dampingRegionOffset;
            simParams->dampingRegionRhs = simParams->drivingRegionRhs + simParams->dampingRegionOffset;
        }
        if ( simParams->numSpinsDampingWidth < 0 ) {
            simParams->numSpinsDampingWidth = simParams->dampingRegionRhs - simParams->dampingRegionLhs + 1;
        }
        simParams->numSpinsDampingGradient = (simParams->numSpinsDampingWidth - simParams->numSpinsDampingPeak) / 2;

        if (simFlags->isDampingLinearAcrossMap) {
            simParams->dampingGradientGradient = simParams->dampingGradientPeak;
        } else {
            simParams->dampingGradientGradient = simParams->gilbertDamping;
        }

        _generateMapOfScaling("Damping Region", simStates->dampingGradientMap, simParams->dampingRegionLhs,
                              simParams->dampingRegionRhs, simParams->numSpinsDampingWidth, simParams->numSpinsDampingGradient,
                              simParams->numSpinsDampingPeak, simParams->dampingGradientGradient, simParams->dampingGradientPeak);
    }
}

void
SolversConfiguration::_generateAbsorbingRegionsUpdated( int numSpinsInChain, int numSpinsAbsorbingRegion,
                                                        double gilbertSpinChain,
                                                        double gilbertAbsorbingRegionInner,
                                                        double gilbertAbsorbingRegionOuter ) {

    // Check for invalid input.
    if ( numSpinsAbsorbingRegion < 0 ) {
        std::cout << "numSpinsAbsorbingRegion is less than zero!";
        exit(1);
    }

    std::vector<double> gilbertMainChain(numSpinsInChain, gilbertSpinChain);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These damping regions are part of my Periodic Boundary Conditions (PBCs) as Absorbing Boundary Conditions (ABCs).
    LinspaceClass DampingRegionLeft;
    LinspaceClass DampingRegionRight;

    DampingRegionLeft.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion,
                                 true, false, false);
    DampingRegionRight.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion,
                                  true, false, false);
    std::vector<double> dampingRegionLHS = DampingRegionLeft.generate_array();
    std::vector<double> dampingRegionRHS = DampingRegionRight.generate_array();

    // std::cout << "dampingRegionLHS.size(): " << dampingRegionLHS.size() << " | dampingRegionRHS.size(): " << dampingRegionRHS.size() << "\n";
    // PrintVector(dampingRegionLHS, false);
    // PrintVector(dampingRegionRHS, false);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Combine all damped regions to form vector which describes the entire spinchain.
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionLHS.begin(), dampingRegionLHS.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), gilbertMainChain.begin(), gilbertMainChain.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionRHS.begin(), dampingRegionRHS.end());
    simStates->gilbertVector.push_back(0);

    // std::cout << "simStates->gilbertVector.size(): " << simStates->gilbertVector.size() << "\n";
    // PrintVector(simStates->gilbertVector, true);
}

void
SolversConfiguration::_generateAbsorbingRegions( int numSpinsInChain, int numSpinsAbsorbingRegion, int numSpinsDRPeak,
                                                 int numSpinsDRGradient, double gilbertSpinChain, double gilbertDRPeak,
                                                 double gilbertAbsorbingRegionInner,
                                                 double gilbertAbsorbingRegionOuter ) {

    // Check for invalid input.
    if ( numSpinsAbsorbingRegion < 0 ) {
        std::cout << "numSpinsAbsorbingRegion is less than zero!";
        exit(1);
    }

    if (numSpinsDRPeak < 0 || numSpinsDRGradient < 0) {
        std::cout << "numSpinsDRPeak or numSpinsDRGradient are less than zero!";
        exit(1);
    }
    int numSpinsChainLhsOfDR = simParams->drivingRegionLhs - numSpinsAbsorbingRegion;
    int numSpinsChainRhsOfDR = numSpinsInChain - simParams->drivingRegionRhs + numSpinsAbsorbingRegion;

    //std::cout << "Before adjustment" << std::endl;
    // std::cout << "numSpinsTotal: " << simParams->systemTotalSpins << " | numSpinsInChain: " << numSpinsInChain << " | numSpinsABC (each): " << simParams->numSpinsInABC << "\n";
    //std::cout << "drivingRegionLhs: " << simParams->drivingRegionLhs << " | drivingRegionRhs: " << simParams->drivingRegionRhs << " | drivingRegionWidth: " << simParams->drivingRegionWidth << "\n";
    //std::cout << "numSpinsChainLhsOfDR: " << numSpinsChainLhsOfDR << " | numSpinsChainRhsOfDR: " << numSpinsChainRhsOfDR << " | numSpinsOutsideOfDR: " << numSpinsChainLhsOfDR + numSpinsChainRhsOfDR << "\n";
    //std::cout << "Spins in chain: " << simParams->numSpinsInChain << " | total sites: " << simParams->systemTotalSpins << " | abc sites: " << simParams->numSpinsInABC  <<std::endl;


    if (simFlags->shouldDriveLHS) {
        // The LHS driving case subtracts 1 from drivingRegionLhs for the zeroth spin addition correction
        numSpinsChainLhsOfDR--; // Just correct for zero-indexing
    } else if (simFlags->shouldDriveRHS) {
        // The RHS driving case subtracts 1 from drivingRegionRhs for the zeroth spin addition correction
        numSpinsChainRhsOfDR--;
    } else if (simFlags->shouldDriveCentre) {
        if ( simParams->drivingRegionWidth % 2 == 0 ) {
            // Even driving region width
            numSpinsChainLhsOfDR--; // Just correct for zero-indexing
        } else {
            // Odd driving region width
            numSpinsChainLhsOfDR--; // We require numRHS > numLHS to maintain some first driven site (Obsidian: 2024-03-01)
        }
    }
        // No adjustment needed for the RHS in the even case, or any case for the odd driving region width
    //std::cout << "after adjustment" << std::endl;
    // std::cout << "numSpinsTotal: " << simParams->systemTotalSpins << " | numSpinsInChain: " << numSpinsInChain << " | numSpinsABC (each): " << simParams->numSpinsInABC << "\n";
    //std::cout << "drivingRegionLhs: " << simParams->drivingRegionLhs << " | drivingRegionRhs: " << simParams->drivingRegionRhs << " | drivingRegionWidth: " << simParams->drivingRegionWidth << "\n";
    //std::cout << "numSpinsChainLhsOfDR: " << numSpinsChainLhsOfDR << " | numSpinsChainRhsOfDR: " << numSpinsChainRhsOfDR << " | numSpinsOutsideOfDR: " << numSpinsChainLhsOfDR + numSpinsChainRhsOfDR << "\n";
    //std::cout << "Spins in chain: " << simParams->numSpinsInChain << " | total sites: " << simParams->systemTotalSpins << std::endl;

    std::vector<double> gilbertMainChain;
    // Calculate the total number of sites in the driving region
    int totalDRWidth = simParams->drivingRegionRhs - simParams->drivingRegionLhs + 1;
    // Calculate expected sizes based on simulation parameters
    int expectedNumGradientSites = 2 * numSpinsDRGradient; // Times two for both left and right gradients
    int expectedNumPeakSites = totalDRWidth - expectedNumGradientSites;

    // Check if the sum of gradient (times 2) and peak sites matches the total driving region width
    if (expectedNumGradientSites + expectedNumPeakSites != totalDRWidth || totalDRWidth != simParams->drivingRegionWidth) {
        throw std::runtime_error("The number of sites in gradients (times 2) and peak does not sum to the total driving region width.");
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Finds the damping regions that are in the main chain, but not part of the driving region width
    LinspaceClass OutsideDrLeft;
    LinspaceClass OutsideDrRight;

    OutsideDrLeft.set_values(gilbertSpinChain, gilbertSpinChain, numSpinsChainLhsOfDR,
                             true, false, false);
    OutsideDrRight.set_values(gilbertSpinChain, gilbertSpinChain, numSpinsChainRhsOfDR,
                              true, false, false);

    std::vector<double> outsideDrLeft = OutsideDrLeft.generate_array();
    std::vector<double> outsideDrRight = OutsideDrRight.generate_array();

    // std::cout << "outsideDrLeft.size(): " << outsideDrLeft.size() << " | outsideDrRight.size(): " << outsideDrRight.size() << "\n";
    // PrintVector(outsideDrLeft, false);
    // PrintVector(outsideDrRight, false);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Finds the damping regions that are in the main chain, and are in the driving region

    std::vector<double> gilbertDROnly;

    LinspaceClass LeftDRGradient;
    LinspaceClass CentreDRPeak;
    LinspaceClass RightDRGradient;

    LeftDRGradient.set_values(gilbertSpinChain, gilbertDRPeak, numSpinsDRGradient,
                              true, false, true);
    CentreDRPeak.set_values(gilbertDRPeak, gilbertDRPeak, numSpinsDRPeak,
                            true, false, false);
    RightDRGradient.set_values(gilbertDRPeak, gilbertSpinChain, numSpinsDRGradient,
                               true, false, true);

    std::vector<double> leftDRGradient = LeftDRGradient.generate_array();
    std::vector<double> centreDRPeak = CentreDRPeak.generate_array();
    std::vector<double> rightDRGradient = RightDRGradient.generate_array();

    // std::cout << "leftDRGradient.size(): " << leftDRGradient.size() << " | centreDRPeak.size(): " << centreDRPeak.size() << " | rightDRGradient.size(): " << rightDRGradient.size() << "\n";
    // PrintVector(leftDRGradient, false);
    // PrintVector(centreDRPeak, false);
    // PrintVector(rightDRGradient, false);

    gilbertDROnly.insert(gilbertDROnly.end(), leftDRGradient.begin(), leftDRGradient.end());
    gilbertDROnly.insert(gilbertDROnly.end(), centreDRPeak.begin(), centreDRPeak.end());
    gilbertDROnly.insert(gilbertDROnly.end(), rightDRGradient.begin(), rightDRGradient.end());

    // std::cout << "gilbertDROnly.size(): " << gilbertDROnly.size() << "\n";
    // PrintVector(gilbertDROnly, false);

    // Check if the lengths of the vectors match the expected gradient/peak widths
    if (leftDRGradient.size() != numSpinsDRGradient || rightDRGradient.size() != numSpinsDRGradient) {
        throw std::runtime_error("The size of gradient alpha vectors does not match the expected gradient width.");
    }

    if (centreDRPeak.size() != expectedNumPeakSites) {
        throw std::runtime_error("The size of the centre peak alpha vector does not match the expected peak width.");
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Combine the main chain all together
    gilbertMainChain.insert(gilbertMainChain.end(), outsideDrLeft.begin(), outsideDrLeft.end());
    gilbertMainChain.insert(gilbertMainChain.end(), gilbertDROnly.begin(), gilbertDROnly.end());
    gilbertMainChain.insert(gilbertMainChain.end(), outsideDrRight.begin(), outsideDrRight.end());

    // std::cout << "gilbertMainChain.size(): " << gilbertMainChain.size() << "\n";
    // PrintVector(gilbertMainChain, false);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Obtain the gradient of the DR region so I can scale other properties (like DMI)
    LinspaceClass LeftDRGradientScaled;
    LinspaceClass CentreDRPeakScaled;
    LinspaceClass RightDRGradientScaled;
    LeftDRGradientScaled.set_values(0.0, 1.0, numSpinsDRGradient, true, false, true);
    CentreDRPeakScaled.set_values(1.0, 1.0, numSpinsDRPeak, true, false, false);
    RightDRGradientScaled.set_values(1.0, 0.0, numSpinsDRGradient, true, false, true);

    std::vector<double> leftDRGradientScaled = LeftDRGradientScaled.generate_array();
    std::vector<double> centreDRPeakScaled = CentreDRPeakScaled.generate_array();
    std::vector<double> rightDRGradientScaled = RightDRGradientScaled.generate_array();
    std::vector<int> leftDRGradientSites;
    std::vector<int> centreDRPeakSites;
    std::vector<int> rightDRGradientSites;

    // std::cout << "leftDRGradientScaled.size(): " << leftDRGradientScaled.size() << " | centreDRPeakScaled: " << centreDRPeakScaled.size() << " | rightDRGradientScaled.size(): " << rightDRGradientScaled.size() << "\n";
    // PrintVector(leftDRGradientScaled, false);
    // PrintVector(centreDRPeakScaled, false);
    // PrintVector(rightDRGradientScaled, false);

    for ( int i = simParams->drivingRegionLhs; i < simParams->drivingRegionLhs + numSpinsDRGradient; i++ ) {
        leftDRGradientSites.push_back(i);
    }
    for ( int i = simParams->drivingRegionLhs + numSpinsDRGradient; i <= simParams->drivingRegionRhs - numSpinsDRGradient; i++ ) {
        centreDRPeakSites.push_back(i);
    }
    for ( int i = simParams->drivingRegionRhs - numSpinsDRGradient + 1; i <= simParams->drivingRegionRhs; i++ ) {
        rightDRGradientSites.push_back(i);
    }

    // Check if the lengths of the vectors match the expected gradient/peak widths
    if (leftDRGradientSites.size() != numSpinsDRGradient || rightDRGradientSites.size() != numSpinsDRGradient) {
        throw std::runtime_error("The size of gradient sites vectors does not match the expected gradient width.");
    }

    if (centreDRPeakSites.size() != expectedNumPeakSites) {
        throw std::runtime_error("The size of the centre peak sites vector does not match the expected peak width.");
    }

    // Check if the lengths of the vectors match the expected gradient/peak widths
    if (leftDRGradientScaled.size() != numSpinsDRGradient || rightDRGradientScaled.size() != numSpinsDRGradient) {
        throw std::runtime_error("The size of lhs/rhs gradient scaled vectors does not match the expected gradient width.");
    }

    if (centreDRPeakScaled.size() != expectedNumPeakSites) {
        throw std::runtime_error("The size of the centre peak scaled vector does not match the expected peak width.");
    }



    // std::cout << "leftDRGradientSites.size(): " << leftDRGradientSites.size() << " | centreDRPeakSites: " << centreDRPeakSites.size() << " | rightDRGradientSites.size(): " << rightDRGradientSites.size() << "\n";
    // PrintVector(leftDRGradientSites, false);
    // PrintVector(centreDRPeakSites, false);
    // PrintVector(rightDRGradientSites, false);

    // Before populating the map, check if sizes match
    if (leftDRGradientSites.size() != leftDRGradientScaled.size()) {
        throw std::runtime_error("Size of leftDRGradientScaled does not match size of leftDRGradientSites");
    }

    // Populate the maps
    for (size_t i = 0; i < leftDRGradientSites.size(); ++i) {
        simStates->dRGradientMap[leftDRGradientSites[i]] = leftDRGradientScaled[i];
    }
    for (size_t i = 0; i < centreDRPeakSites.size(); ++i) {
        simStates->dRGradientMap[centreDRPeakSites[i]] = centreDRPeakScaled[i];
    }
    for (size_t i = 0; i < rightDRGradientSites.size(); ++i) {
        simStates->dRGradientMap[rightDRGradientSites[i]] = rightDRGradientScaled[i];
    }
    // std::cout << "simStates->dRGradientMap.size(): " << simStates->dRGradientMap.size() << "\n";
    // PrintVector(simStates->dRGradientMap, false);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These damping regions are part of my Periodic Boundary Conditions (PBCs) as Absorbing Boundary Conditions (ABCs).
    LinspaceClass DampingRegionLeft;
    LinspaceClass DampingRegionRight;

    DampingRegionLeft.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion,
                                 true, false, false);
    DampingRegionRight.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion,
                                  true, false, false);
    std::vector<double> dampingRegionLHS = DampingRegionLeft.generate_array();
    std::vector<double> dampingRegionRHS = DampingRegionRight.generate_array();

    // std::cout << "dampingRegionLHS.size(): " << dampingRegionLHS.size() << " | dampingRegionRHS.size(): " << dampingRegionRHS.size() << "\n";
    // PrintVector(dampingRegionLHS, false);
    // PrintVector(dampingRegionRHS, false);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Combine all damped regions to form vector which describes the entire spinchain.
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionLHS.begin(), dampingRegionLHS.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), gilbertMainChain.begin(), gilbertMainChain.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionRHS.begin(), dampingRegionRHS.end());
    simStates->gilbertVector.push_back(0);

    // std::cout << "simStates->gilbertVector.size(): " << simStates->gilbertVector.size() << "\n";
    // PrintVector(simStates->gilbertVector, true);
}

void SolversConfiguration::_generateMapOfScaling( const std::string &regionName,
                                                  std::map<int, double> &assignedMap,
                                                  const int &regionLhsSite, const int &regionRhsSite,
                                                  const int &regionNumSites, const int &numSitesGradientPartition,
                                                  const int &numSitesPeakPartition, const double &scalingOuter,
                                                  const double &scalingInner ) {

    if (regionLhsSite < 0 || regionRhsSite < 0 || regionNumSites < 0 || numSitesGradientPartition < 0 || numSitesPeakPartition < 0) {
        throw std::runtime_error("[SolversConfiguration>_generateMapOfScaling]: One of the input parameters is less than "
                                 "zero while generating a map for " + regionName);
    }

    if (numSitesGradientPartition * 2 + numSitesPeakPartition != regionNumSites) {
        throw std::runtime_error("[SolversConfiguration>_generateMapOfScaling]: The sum of the gradient and peak sites "
                                 "does not match the total region width of " + regionName + " | "
                                 + std::to_string(numSitesGradientPartition) + " | "
                                 + std::to_string(numSitesPeakPartition) + " | "
                                 + std::to_string(regionNumSites));
    }

    // Generates the map of the DMI vector
    LinspaceClass LeftGradientScaled;
    LinspaceClass CentreScaled;
    LinspaceClass RightGradientScaled;

    LeftGradientScaled.set_values(scalingOuter, scalingInner, numSitesGradientPartition, true, false, true);
    CentreScaled.set_values(scalingInner, scalingInner, numSitesPeakPartition, true, false, false);
    RightGradientScaled.set_values(scalingInner, scalingOuter, numSitesGradientPartition, true, false, true);

    std::vector<double> leftGradientScaled = LeftGradientScaled.generate_array();
    std::vector<double> centreScaled = CentreScaled.generate_array();
    std::vector<double> rightGradientScaled = RightGradientScaled.generate_array();
    std::vector<int> leftGradientSiteIndexes;
    std::vector<int> centreSiteIndexes;
    std::vector<int> rightGradientSiteIndexes;

    // std::cout << "leftDRGradientScaled.size(): " << leftGradientScaled.size() << " | centreDRPeakScaled: " << centreScaled.size() << " | rightDRGradientScaled.size(): " << rightGradientScaled.size() << "\n";
    // PrintVector(leftGradientScaled, false);
    // PrintVector(centreScaled, false);
    // PrintVector(rightGradientScaled, false);

    int leftGradientLhs = regionLhsSite, leftGradientRhs = regionLhsSite + numSitesGradientPartition;
    int centreLhs = leftGradientRhs, centreRhs = regionRhsSite - numSitesGradientPartition;
    int rightGradientLhs = centreRhs + 1, rightGradientRhs = regionRhsSite;

    for ( int i = leftGradientLhs; i < leftGradientRhs; i++ ) {
        leftGradientSiteIndexes.push_back(i);
    }
    for ( int i = centreLhs; i <= centreRhs; i++ ) {
        centreSiteIndexes.push_back(i);
    }
    for ( int i = rightGradientLhs; i <= rightGradientRhs; i++ ) {
        rightGradientSiteIndexes.push_back(i);
    }

    // Check if the lengths of the vectors match the expected gradient/peak widths
    if (leftGradientSiteIndexes.size() != numSitesGradientPartition || rightGradientSiteIndexes.size() != numSitesGradientPartition) {
        throw std::runtime_error("The size of gradient sites vectors does not match the expected gradient width.");
    }

    if ( centreSiteIndexes.size() != numSitesPeakPartition) {
        throw std::runtime_error("The size of the centre peak sites vector does not match the expected peak width.");
    }

    // std::cout << "leftGradientSiteIndexes: " << leftGradientSiteIndexes.size() << " | centreDRPeakSites: " << centreSiteIndexes.size() << " | rightGradientSiteIndexes: " << rightGradientSiteIndexes.size() << "\n";
    // PrintVector(leftGradientSiteIndexes, false);
    // PrintVector(centreSiteIndexes, false);
    // PrintVector(rightGradientSiteIndexes, false);

    // Before populating the map, check if sizes match
    if (leftGradientSiteIndexes.size() != leftGradientScaled.size())
        throw std::runtime_error("Size of leftGradientSiteIndexes does not match size of leftGradientScaled");
    if (centreSiteIndexes.size() != centreScaled.size())
        throw std::runtime_error("Size of centreSiteIndexes does not match size of centreScaled");
    if (rightGradientSiteIndexes.size() != rightGradientScaled.size())
        throw std::runtime_error("Size of rightGradientSiteIndexes does not match size of rightGradientScaled");


    // Populate the maps
    for (size_t i = 0; i < leftGradientSiteIndexes.size(); i++) {
        assignedMap[leftGradientSiteIndexes[i]] = leftGradientScaled[i];
    }
    for (size_t i = 0; i < centreSiteIndexes.size(); i++) {
        assignedMap[centreSiteIndexes[i]] = centreScaled[i];
    }
    for (size_t i = 0; i < rightGradientSiteIndexes.size(); i++) {
        assignedMap[rightGradientSiteIndexes[i]] = rightGradientScaled[i];
    }

    // std::cout << "assignedMap.size(): " << assignedMap.size() << "\n";
    // PrintVector(assignedMap, false);
}
std::vector<double> SolversConfiguration::_generateDampingForSubset( int numSitesInSubset, int numSitesInSubsetCentre,
                                                                     int numSpinsInSubsetLHS, int numSpinsInSubsetRHS,
                                                                     double gilbertSubsetCentre,
                                                                     double gilbertSubsetLHS, double gilbertSubsetRHS,
                                                                     bool includeLeftEndpointZero,
                                                                     bool includeRightEndpointZero ) {
    if (numSitesInSubset != numSitesInSubsetCentre + numSpinsInSubsetLHS + numSpinsInSubsetRHS) {
        std::cout << "numSitesInSubset does not match the sum of the other parameters!";
        exit(1);
    }

    if (numSitesInSubset < 0 || numSitesInSubsetCentre < 0 || numSpinsInSubsetLHS < 0 || numSpinsInSubsetRHS < 0) {
        std::cout << "Invalid input for numSitesInSubset, numSitesInSubsetCentre, numSpinsInSubsetLHS or "
                     "numSpinsInSubsetRHS!";
        exit(1);
    }

    std::vector<double> subsetRegion;

    LinspaceClass CentreRegionLeft;
    LinspaceClass CentreRegionCentre;
    LinspaceClass CentreRegionRight;

    CentreRegionLeft.set_values(gilbertSubsetLHS, gilbertSubsetCentre, numSpinsInSubsetLHS,
                                true, false, false);
    CentreRegionCentre.set_values(gilbertSubsetCentre, gilbertSubsetCentre, numSitesInSubsetCentre,
                                  true, false, false);
    CentreRegionRight.set_values(gilbertSubsetCentre, gilbertSubsetRHS, numSpinsInSubsetRHS,
                                 true, false, false);

    std::vector<double> centreRegionLeft = CentreRegionLeft.generate_array();
    std::vector<double> centreRegionCentre = CentreRegionCentre.generate_array();
    std::vector<double> centreRegionRight = CentreRegionRight.generate_array();

    if (includeLeftEndpointZero)
        subsetRegion.push_back(0);

    subsetRegion.insert(subsetRegion.end(), centreRegionLeft.begin(), centreRegionLeft.end());
    subsetRegion.insert(subsetRegion.end(), centreRegionCentre.begin(), centreRegionCentre.end());
    subsetRegion.insert(subsetRegion.end(), centreRegionRight.begin(), centreRegionRight.end());

    if (includeRightEndpointZero)
        subsetRegion.push_back(0);

    return subsetRegion;
}


void SolversConfiguration::_generateAbsorbingRegions( int numSpinsInChain, int numSpinsAbsorbingRegion,
                                                      double gilbertSpinChain,
                                                      double gilbertAbsorbingRegionInner,
                                                      double gilbertAbsorbingRegionOuter ) {
    // Generate the damping regions that are appended to either end of the spin chain.

    // These damping regions are part of my Periodic Boundary Conditions (PBCs) as Absorbing Boundary Conditions (ABCs).
    LinspaceClass DampingRegionLeft;
    LinspaceClass DampingRegionRight;

    if ( numSpinsAbsorbingRegion < 0 ) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(1);
    }

    std::vector<double> gilbertChain(numSpinsInChain, gilbertSpinChain);

    DampingRegionLeft.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion,
                                 true, false, false);
    DampingRegionRight.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion,
                                  true, false, false);
    std::vector<double> dampingRegionLHS = DampingRegionLeft.generate_array();
    std::vector<double> dampingRegionRHS = DampingRegionRight.generate_array();

    // Combine all damped regions to form vector which describes the entire spinchain.
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionLHS.begin(), dampingRegionLHS.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());
    simStates->gilbertVector.insert(simStates->gilbertVector.end(), dampingRegionRHS.begin(), dampingRegionRHS.end());
    simStates->gilbertVector.push_back(0);
}

void
SolversConfiguration::_setupDrivingRegion( int numSpinsInChain, int numSpinsAbsorbingRegion, int drivingRegionWidth ) {

    if ( simFlags->shouldDriveDiscreteSites ) {
        for ( int &drivenSite: simStates->discreteDrivenSites )
            drivenSite -= 1;

        // Update so that output file is intuitive
        simParams->drivingRegionWidth = 0;
        simParams->drivingRegionLhs = 0;
        simParams->drivingRegionRhs = 0;
        return;
    }

    if ( simFlags->shouldDriveCentre ) {
        // Use same midpoint for even and odd lengthened chains (numSpinsInChain)
        /*
         * CONVENTION:
         *      - The driving region is centred around the midpoint of `floor(numSpinsInChain / 2).
         *      - We want (N) and (N+1) cases to have the same left-hand site in the driving region for consistency,
         *        where N is the number of spins in the chain (numSpinsInChain), and N >= 1.
         *
         *      WHEN numSpinsInChain is EVEN:
         *      - If the driving region width is even, then there are an EQUAL number of sites on the LHS and RHS
         *        of the driving region itself. We ADD 1 to the LHS to ensure that the driving region is centred and
         *        indexed correctly.
         *      - If the driving region width is odd, then there are an UNEQUAL number of sites on the LHS and RHS; one
         *        more site on the RHS than the LHS. We therefore do not need to add 1 to the LHS to adjust for indexing,
         *        as the midpoint is already centred.
         *
         *      WHEN numSpinsInChain is ODD:
         *      - First, to ensure that the first driven site is consistent with the even case, we ADD 1 to the LHS and
         *        to the RHS. This will ensure that indexing adjustments are consistent with the other changes.
         *      - If the driving region width is even, then there are an UNEQUAL number of sites on the LHS and RHS; one
         *        more site on the LHS than the RHS. We therefore ADD 1 to the LHS to ensure that the LHS has more sites
         *        than the RHS, and that the driving region is centred and indexed correctly.
         *      - If the driving region width is odd, then there are an EQUAL number of sites on the LHS and RHS of the
         *        driving region itself. We do not need to add 1 to the LHS to adjust for indexing, as the midpoint is
         *        already centred.
         */
        if ( drivingRegionWidth % 2 == 0 ) {
            // By convention even case shortens the LHS
            simParams->drivingRegionLhs =
                    (numSpinsInChain / 2) + numSpinsAbsorbingRegion - (drivingRegionWidth / 2) + 1;
            simParams->drivingRegionRhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion + (drivingRegionWidth / 2);
        } else {
            // Use convention to round down then add additional site (odd number) to the RHS of the centre
            simParams->drivingRegionLhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion - (drivingRegionWidth / 2);
            simParams->drivingRegionRhs = (numSpinsInChain / 2) + numSpinsAbsorbingRegion + (drivingRegionWidth / 2);
        }

        if ( simParams->numSpinsInChain % 2 == 1) {
            simParams->drivingRegionLhs++;
            simParams->drivingRegionRhs++;
        }
        return;
    }

    if ( simFlags->shouldDriveBothSides ) {
        throw std::runtime_error(
                std::string("shouldDriveBothSides has been selected, but this feature is not yet implemented"));
    }

    if ( simFlags->shouldDriveLHS ) {
        // The +1/-1 offset excludes the zeroth spin while retaining the correct driving width
        simParams->drivingRegionLhs = numSpinsAbsorbingRegion + 1;
        simParams->drivingRegionRhs = simParams->drivingRegionLhs + drivingRegionWidth - 1;
        return;
    }

    if ( simFlags->shouldDriveRHS ) {
        // The +1 is to correct the offset of adding a zeroth spin
        simParams->drivingRegionRhs = simParams->systemTotalSpins - numSpinsAbsorbingRegion - 1;
        simParams->drivingRegionLhs = simParams->drivingRegionRhs - drivingRegionWidth + 1;
        return;
    }

    if ( simFlags->hasCustomDrivePosition ) {
        // The +1/-1 offset excludes the zeroth spin while retaining the correct driving width
        simParams->drivingRegionLhs = simParams->numSpinsInABC + 1;
        simParams->drivingRegionRhs = simParams->drivingRegionLhs + drivingRegionWidth - 1;
        return;
    }

}

void SolversConfiguration::_generateExchangeVector( int numSpinsAbsorbingRegion, int numSpinPairs, double exchangeMin,
                                                    double exchangeMax ) {
    /*
     * Create the arrays which house the exchange integral values. There are options to have a non-uniform exchange coded in, as well as the option to
     * induce a 'kick' into the system by initialising certain spins to have differing parameters to their neighbours.
     */
    LinspaceClass SpinChainExchange;
    LinspaceClass SpinChainExchangeLHS;
    LinspaceClass SpinChainExchangeRHS;

    if ( numSpinsAbsorbingRegion > 0 ) {
        if ( simFlags->hasSingleExchangeRegion ) {
            numSpinPairs += (2 * simParams->numSpinsInABC);
            SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, true, false);
            simStates->exchangeVec = SpinChainExchange.generate_array();
        } else {
            simStates->exchangeVec.push_back(0.0);

            SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, false, false);
            std::vector<double> tempExchangeChain = SpinChainExchange.generate_array();

            // SpinChainExchangeLHS.set_values(exchangeMin, exchangeMax, numSpinsAbsorbingRegion, true, false);
            // std::vector<double> tempExchangeLHS = SpinChainExchangeLHS.generate_array();
            std::vector<double> tempExchangeLHS(numSpinsAbsorbingRegion, exchangeMin);

            // SpinChainExchangeRHS.set_values(exchangeMax, exchangeMin, numSpinsAbsorbingRegion, true, false);
            // std::vector<double> tempExchangeRHS = SpinChainExchangeRHS.generate_array();
            std::vector<double> tempExchangeRHS(numSpinsAbsorbingRegion, exchangeMax);

            simStates->exchangeVec.insert(simStates->exchangeVec.end(), tempExchangeLHS.begin(), tempExchangeLHS.end());
            simStates->exchangeVec.insert(simStates->exchangeVec.end(), tempExchangeChain.begin(),
                                          tempExchangeChain.end());
            simStates->exchangeVec.insert(simStates->exchangeVec.end(), tempExchangeRHS.begin(), tempExchangeRHS.end());
            simStates->exchangeVec.push_back(0.0);
        }
    } else {
        // The linearly spaced vector is saved as the class member 'exchangeVec' simply to increase code readability
        SpinChainExchange.set_values(exchangeMin, exchangeMax, numSpinPairs, true, true, false);
        simStates->exchangeVec = SpinChainExchange.generate_array();
    }
}

void SolversConfiguration::_setupInitMagneticMoments( double mxInit, double myInit, double mzInit ) {

    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mxInitCond(simParams->systemTotalSpins, mxInit), myInitCond(simParams->systemTotalSpins,
                                                                                    myInit), mzInitCond(
            simParams->systemTotalSpins, mzInit);
    // mxInitCond[0] = mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < systemTotalSpins; i++) {
        InitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    if ( !simFlags->isFerromagnetic ) {
        for ( int i = 0; i < simParams->systemTotalSpins; i++ ) {
            if ( i % 2 == 1 )
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

void SolversConfiguration::resetInitMagneticMoments( double mxInit, double myInit, double mzInit ) {

    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions

    std::fill(simStates->mx0.begin(), simStates->mx0.end(), mxInit);
    std::fill(simStates->my0.begin(), simStates->my0.end(), myInit);
    std::fill(simStates->mz0.begin(), simStates->mz0.end(), mzInit);

    if ( !simStates->mx0.empty()) {
        int vecLength = simStates->mx0.size() - 1;
        simStates->mx0[0] = 0;
        simStates->mx0[vecLength] = 0;

        simStates->my0[0] = 0;
        simStates->my0[vecLength] = 0;

        simStates->mz0[0] = 0;
        simStates->mz0[vecLength] = 0;
    }

    // mxInitCond[0] = mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < systemTotalSpins; i++) {
        InitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    if ( !simFlags->isFerromagnetic ) {
        for ( int i = 0; i < simParams->systemTotalSpins; i++ ) {
            if ( i % 2 == 1 )
                simStates->mz0[i] *= -1.0;
        }
    }

}

void SolversConfiguration::_testShockwaveInitConditions() {

    if ( simFlags->hasShockwave ) {
        simParams->shockwaveStepsize =
                (simParams->shockwaveMax - simParams->shockwaveInitialStrength) / simParams->shockwaveGradientTime;
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

void SolversConfiguration::_generateMultilayerAbsorbingRegions( int numSpinsAbsorbingRegion, double gilbertSpinChain,
                                                                double gilbertAbsorbingRegionInner,
                                                                double gilbertAbsorbingRegionOuter ) {

    LinspaceClass AbsorbingRegionLHS;
    LinspaceClass AbsorbingRegionRHS;

    if ( numSpinsAbsorbingRegion < 0 ) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(1);
    }

    for ( int &layer: simStates->layerSpinsInChain ) {
        std::vector<double> gilbertChain(layer, gilbertSpinChain);

        AbsorbingRegionLHS.set_values(gilbertAbsorbingRegionOuter, gilbertAbsorbingRegionInner, numSpinsAbsorbingRegion,
                                      true, false, false);
        AbsorbingRegionRHS.set_values(gilbertAbsorbingRegionInner, gilbertAbsorbingRegionOuter, numSpinsAbsorbingRegion,
                                      true, false, false);

        // Change of name to be specific that this is only the damping region
        std::vector<double> dampingRegionLHS = AbsorbingRegionLHS.generate_array();
        std::vector<double> dampingRegionRHS = AbsorbingRegionRHS.generate_array();

        // Combine all damped regions to form vector which describes the entire spinchain.
        simStates->gilbertVectorMulti[layer].insert(simStates->gilbertVectorMulti[layer].end(),
                                                    dampingRegionLHS.begin(), dampingRegionLHS.end());
        simStates->gilbertVectorMulti[layer].insert(simStates->gilbertVectorMulti[layer].end(), gilbertChain.begin(),
                                                    gilbertChain.end());
        simStates->gilbertVectorMulti[layer].insert(simStates->gilbertVectorMulti[layer].end(),
                                                    dampingRegionRHS.begin(), dampingRegionRHS.end());
        simStates->gilbertVectorMulti[layer].push_back(0);

        //PrintVector(gilbertVectorMulti[i], false);
    }
}

void SolversConfiguration::_setupInitMultilayerMagneticMoments(
        std::vector<std::vector<std::vector<double>>> &nestedNestedVector,
        int layer, double mxInit, double myInit, double mzInit ) {

    // mxInitCond[0] = mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < systemTotalSpins; i++) {
        mxInitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    for ( int i = 0; i < simStates->layerTotalSpins[layer]; i++ ) {
        // i is simply an index here
        nestedNestedVector[layer].push_back({mxInit, myInit, mzInit});
    }

    // This zero is the (N+1)th spin on the RHS of the chain
    nestedNestedVector[layer].push_back({0.0, 0.0, 0.0});

}

std::vector<std::vector<std::vector<double>>> SolversConfiguration::initializeNestedNestedVector( int numSites,
                                                                                                  bool includeEnd ) {
    /* Legacy code not used currently. Example implementation is below
     * std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested3;
     * mValsNested3["nestedNestedVector3"] = initializeNestedNestedVector(1, true);
     * std::vector<std::vector<std::vector<double>>> m2Nest = mValsNested3["nestedNestedVector3"];
     * _setupInitMultilayerMagneticMoments(m2Nest, 1, 0, 0 , 0);
     */
    std::vector<std::vector<std::vector<double>>> innerNestedVector;
    for ( int j = 0; j < numSites; j++ ) {
        std::vector<std::vector<double>> innerVector;
        std::vector<double> innermostVector = {0.0, 0.0, 0.0};
        innerVector.push_back(innermostVector);
        innerNestedVector.push_back(innerVector);
    }
    return innerNestedVector;
}

std::vector<std::vector<std::vector<double>>>
SolversConfiguration::InitialiseNestedVectors( int &totalLayer, double &mxInit, double &myInit, double &mzInit ) {

    // Initialise mapping
    std::map<std::string, std::vector<std::vector<std::vector<double>>>> mTermsMapping;

    // This is likely a very slow way to initialise (push_back is slow), but this works for now. Fix if it is a bottleneck later
    std::vector<std::vector<std::vector<double>>> innerNestedVector;
    for ( int i = 0; i < totalLayer; i++ ) {
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
    for ( int layer = 0; layer < totalLayer; layer++ )
        _setupInitMultilayerMagneticMoments(mTermsNested, layer, mxInit, myInit, mzInit);

    return mTermsNested;
}