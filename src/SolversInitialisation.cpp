//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/SolversInitialisation.h"

SolversInitialisation::SolversInitialisation(std::shared_ptr<SimulationManager> sharedSimManager,
                                             std::shared_ptr<SimulationParameters> sharedSimParams,
                                             std::shared_ptr<SimulationStates> sharedSimStates,
                                             std::shared_ptr<SimulationFlags> sharedSimFlags)

    : SolversSuperClass(std::move(sharedSimManager), std::move(sharedSimParams), std::move(sharedSimStates),  std::move(sharedSimFlags)) {
}

void SolversInitialisation::Initialise() {
    // TODO. Add method that automatically calls components which contain preset values (e.g. for maced2022breaking Fig.2)
    // Keep this set order to avoid errors
    _setSimulationFlags();
    _setSimulationParameters();
    _generateRemainingParameters();
    _setMaterialParameters();
    _guardClauses();
}
void SolversInitialisation::_recalculateParameters() {
    // Keep this set order to avoid errors
    _generateRemainingParameters();
    _setMaterialParameters();
    _guardClauses();
}

void SolversInitialisation::_setSimulationFlags() {

    // Debugging Flags
    simFlags->shouldTrackMagneticMomentNorm = true;

    // Model Type
    simFlags->shouldUseLLG = true;
    simFlags->shouldUseSLLG = false;

    // Magnetic Interaction Flags
    simFlags->hasRisingTime = false;
    simFlags->hasDipolar = false;
    simFlags->hasDMI = true;
    simFlags->hasSTT = false;
    simFlags->hasStaticZeeman = true;
    simFlags->hasDmi1DThinFilm = false;
    simFlags->hasDemagFactors = false;
    simFlags->hasDemagFFT = false;
    simFlags->hasShapeAnisotropy = false;

    // Material Flags
    simFlags->hasMultipleAppendedLayers = false;
    simFlags->hasMultipleStackedLayers = false;
    simFlags->hasSingleExchangeRegion = true;

    // Drive Position Flags
    simFlags->shouldDriveDiscreteSites = false;
    simFlags->hasCustomDrivePosition = false;
    simFlags->shouldDriveAllLayers = true;
    simFlags->shouldDriveBothSides = false;
    simFlags->shouldDriveCentre = true;
    simFlags->shouldDriveLHS = false;
    simFlags->shouldDriveRHS = false;

    // Drive Manipulation Flags
    simFlags->isOscillatingZeemanStatic = false;
    simFlags->shouldDriveCeaseEarly = false;

    // Output Flags
    simFlags->shouldPrintAllData = false;
    simFlags->shouldPrintDiscreteTimes = true;
    simFlags->shouldPrintDiscreteSites = false;

    // TESTING ONLY
    simFlags->hasGradientRegionForOscillatingZeeman = false;
    simFlags->hasGradientRegionForDmi = false;
    simFlags->hasGradientRegionForDamping = false;

    simFlags->shouldDmiGradientMirrorOscillatingZeeman = false;
    simFlags->shouldDampingGradientMirrorOscillatingZeeman = false;

    simFlags->shouldRestrictDmiToWithinGradientRegion = false;
    simFlags->isOscillatingZeemanLinearAcrossMap = false;
    simFlags->isDmiLinearAcrossMap = false;
    simFlags->isDampingLinearAcrossMap = false;

    simFlags->useGenerateABCUpdated = true;
}

void SolversInitialisation::_setSimulationParameters() {
    // Old global variables that have now been moved here
    simParams->exchangeInteraction.setUniformExchange(40, CommonStructures::Unit::Tesla);
    // exchangeInteraction.setNonUniformExchange(132.5, 132.5, CommonStructures::Unit::Tesla);
    simParams->exchangeEnergyMin = simParams->exchangeInteraction.uniformStrength;
    simParams->exchangeEnergyMax = simParams->exchangeInteraction.uniformStrength;

    // THIS MUST BE MANUALLY SET FOR NOW
    simParams->latticeConstant = 1e-9;

    simFlags->equilibriumOrientation = Axis::y;

    // Main Parameters
    simParams->ambientTemperature = 0; // Kelvin
    simParams->drivingFreq = 20.0 * 1e9;
    simParams->staticZeemanStrength[simFlags->equilibriumOrientation] = 1e-2;
    simParams->oscillatingZeemanStrength = 1e-2;
    simParams->forceStopAtIteration = -1;
    simParams->gyroMagConst = GV.GetGyromagneticConstant();
    simParams->maxSimTime = 1e-9;
    simParams->satMag = 800e3;
    simParams->stepsize = 1e-15;

    // Data Output Parameters
    simStates->fixedOutputSites = {1500};
    simParams->numberOfDataPoints = 1e2; //static_cast<int>(maxSimTime / recordingInterval);

    // Damping Factors
    simParams->gilbertDamping.factor = 1e-4;
    simParams->gilbertDamping.innerABC = 1e-4;
    simParams->gilbertDamping.outerABC = 1e0;

    // Spin chain and multi-layer Parameters
    simStates->discreteDrivenSites = {1};
    simParams->drivingRegion.regionWidth = 200;
    simParams->drivingRegion.peakWidth = simParams->drivingRegion.regionWidth;

    simParams->numNeighbours = -1;
    simParams->numSpinsInABC = 300;
    simParams->numLayers = 1;

    // 'Rising time' parameters
    simParams->risingTimeStartAtIteration = 0.0;
    simParams->risingTimePeriod = 1;
    simParams->risingTimeInitialMagnitude = 0;  // Set equal to oscillatingZeemanStrength if NOT starting at time=0
    simParams->risingTimeMaximum = simParams->oscillatingZeemanStrength;
    simParams->risingTimeScalingFactor = 1;

    // Currently unused parameters
    simParams->recordingInterval = 2.5e-15;
    simParams->layerOfInterest = 1;

    // TODO Transferred from _generateRemainingParameters for now

    // Set flags from global variables
    simFlags->isFerromagnetic = GV.GetIsFerromagnetic();
    simParams->dmiConstant = GV.GetDMIConstant() / 2; // Needed to match my Python and C++ code
    simParams->dmiVector.x() = simParams->dmiConstant;

    simParams->anisotropyField = GV.GetAnisotropyField();
    simParams->numSpinsInChain = GV.GetNumSpins();
    simParams->exchangeStiffness = 5.3e-17;  // TODO. Integrate with ExchangeMinVal parameter
    simParams->spinPolarisation = 0.5;
    simParams->spinTransferEfficiency = 0.4;
    simParams->PERMEABILITY_IRON *= _BOHR_MAGNETON;  // Conversion to Am^2
    simParams->dipoleConstant = SimulationParameters::PERM_FREESPACE / (4.0 * M_PI);
    simParams->layerOfInterest -= 1;  // To correct for 0-indexing

    // Testing ONLY!
    simParams->drivingRegion.gradientWidth =  (simParams->drivingRegion.regionWidth - simParams->drivingRegion.peakWidth) / 2;

    simParams->dmiRegion.offsetWidth= 300;
    simParams->dmiRegion.peakWidth = simParams->drivingRegion.regionWidth;

    simParams->dampingRegion.offsetWidth = simParams->dmiRegion.offsetWidth;
    simParams->dampingRegion.peakValue = 1e-2;
    simParams->dampingRegion.peakWidth = simParams->drivingRegion.regionWidth;
}

void SolversInitialisation::_generateRemainingParameters() {

    // Computations based upon other inputs
    simParams->drivingAngFreq = 2 * M_PI * simParams->drivingFreq;

    simParams->risingTimePeriod /= simParams->stepsize;
    simParams->risingTimeEndAtIteration = (simParams->risingTimePeriod / simParams->maxSimTime);

    simParams->iterationEnd = static_cast<int>(simParams->maxSimTime / simParams->stepsize);
    simParams->stepsizeHalf = simParams->stepsize / 2.0;

    simStates->layerSpinsInChain = {simParams->drivingRegion.regionWidth, simParams->numSpinsInChain};
    simStates->layerSpinPairs.clear();
    simStates->layerTotalSpins.clear();
    for (int& spinsInChain: simStates->layerSpinsInChain) {
        simStates->layerSpinPairs.push_back(spinsInChain - 1);
        simStates->layerTotalSpins.push_back(spinsInChain + 2 * simParams->numSpinsInABC);
    }

    simStates->gilbertVectorMulti.resize(simParams->numLayers, {0});
    simParams->numberOfSpinPairs = simParams->numSpinsInChain - 1;

    simParams->systemTotalSpins = simParams->numSpinsInChain + 2 * simParams->numSpinsInABC;
}

void SolversInitialisation::_setMaterialParameters() {

    if (simFlags->isFerromagnetic)
        simParams->anisotropyField = 0;

    if (!simFlags->hasStaticZeeman)
        simParams->staticZeemanStrength.z() = 0.0;
}

void SolversInitialisation::_guardClauses() {

    if (simFlags->shouldDriveCeaseEarly and simParams->risingTimeEndAtIteration <= 0) {
        std::cout << "Warning: [shouldDriveCeaseEarly: True] however [risingTimeEndAtIteration: " << simParams->risingTimeEndAtIteration << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (simFlags->hasRisingTime and simParams->risingTimeStartAtIteration < 0) {
        std::cout << "Warning: [hasRisingTime: True] however [risingTimeStartAtIteration: " << simParams->risingTimeStartAtIteration << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if ((simFlags->shouldPrintDiscreteSites and simFlags->shouldPrintDiscreteTimes) or (simFlags->shouldPrintDiscreteSites and simFlags->shouldPrintAllData) or
        (simFlags->shouldPrintDiscreteTimes and simFlags->shouldPrintAllData)) {
        std::cout << "Warning: Multiple output flags detected. [shouldPrintDiscreteSites: " << simFlags->shouldPrintDiscreteSites
                  << "] | [shouldPrintDiscreteTimes: " << simFlags->shouldPrintDiscreteTimes << "] | [shouldPrintAllData: " << simFlags->shouldPrintAllData << "]"
                  << std::endl;
        exit(1);
    }

    if ((simFlags->shouldDriveLHS && simFlags->shouldDriveCentre) || (simFlags->shouldDriveLHS && simFlags->shouldDriveBothSides)
        || (simFlags->shouldDriveCentre && simFlags->shouldDriveBothSides) || (simFlags->shouldDriveLHS && simFlags->shouldDriveRHS)
        || (simFlags->shouldDriveCentre && simFlags->shouldDriveRHS) || (simFlags->shouldDriveRHS && simFlags->shouldDriveBothSides)) {
        std::cout << "Warning: two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << simFlags->shouldDriveLHS << "\n_centralDrive: " << simFlags->shouldDriveCentre
                  << "\n_dualDrive: " << simFlags->shouldDriveBothSides << "\n _rhsDrive: " << simFlags->shouldDriveRHS << "\n\nExiting...";
        exit(1);
    }

    if (simFlags->shouldPrintDiscreteSites and simStates->fixedOutputSites.empty()) {
        std::cout << "Warning: Request to print fixed sites, but no sites were given [fixedOutputSites: (";
        for (int & fixed_out_val : simStates->fixedOutputSites)
                std::cout << fixed_out_val << ", ";
        std::cout << ")].";
        exit(1);
    }

    if (simParams->numberOfDataPoints > simParams->iterationEnd) {
        std::cout << "Warning: You tried to print more data than was generated [numberOfDataPoints > iterationEnd]";
        exit(1);
    }

    if (simFlags->shouldUseLLG and simFlags->shouldUseSLLG) {
        std::cout << "Warning: You cannot use both the magDynamics and sLLG equations. Please choose one or the other.";
        exit(1);
    }

    if ( simFlags->hasMultipleStackedLayers and simParams->numLayers < 2) {
        std::cout << "Warning: You cannot use the multilayer solver with less than 2 layers.";
        exit(1);
    }

    if (simFlags->hasDemagFactors && simFlags->hasDemagFFT) {
        std::cout << "Warning: You cannot use both the intense and FFT demag solvers. Please choose one or the other.";
        exit(1);
    }

    if ((simFlags->hasDemagFactors && !simFlags->isFerromagnetic) || (simFlags->hasDemagFFT && !simFlags->isFerromagnetic)) {
        std::cout << "Warning: You cannot use the demag solvers with non-ferromagnetic materials.";
        exit(1);
    }

    if (simFlags->shouldDriveDiscreteSites && simStates->discreteDrivenSites.empty()) {
        std::cout << "Warning: Request to drive discrete sites, but no sites were given [discreteDrivenSites].";
        exit(1);
    }

    if (simFlags->hasShapeAnisotropy && simFlags->isFerromagnetic) {
        std::cout << "Warning: You cannot use shape anisotropy with ferromagnetic materials.";
        exit(1);
    }

    if ( (simParams->drivingRegion.peakWidth + simParams->drivingRegion.gradientWidth) > simParams->drivingRegion.regionWidth) {
        std::cout << "Warning: Attempted DR gradient setup exceeds the drivingRegionWidth!";
        exit(1);
    }

}

void SolversInitialisation::testModifyingDouble(double newValue) {
    // Test modifying a double in the parent class from the child
    simParams->ambientTemperature = newValue;
    simParams->risingTimeEndAtIteration = 0.12345;
    std::cout << "I changed another value: " << simParams->risingTimeEndAtIteration << std::endl;
}