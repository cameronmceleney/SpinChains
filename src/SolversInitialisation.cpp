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
    simFlags->hasShockwave = false;
    simFlags->hasDipolar = false;
    simFlags->hasDMI = true;
    simFlags->hasSTT = false;
    simFlags->hasStaticZeeman = true;
    simFlags->hasDemag1DThinFilm = false;
    simFlags->hasDemagIntense = false;
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
    simFlags->hasGradientRegionForOscillatingZeeman = true;
    simFlags->hasGradientRegionForDmi = true;
    simFlags->hasGradientRegionForDamping = true;

    simFlags->shouldDmiGradientMirrorOscillatingZeeman = true;
    simFlags->shouldDampingGradientMirrorOscillatingZeeman = true;

    simFlags->shouldRestrictDmiToWithinGradientRegion = true;
    simFlags->hasOscillatingZeemanSharpInterface = true;

    simFlags->useGenerateABCUpdated = true;
}

void SolversInitialisation::_setSimulationParameters() {

    simParams->latticeConstant = 1e-9;  // THIS MUST BE MANUALLY SET FOR NOW
    // Main Parameters
    simParams->ambientTemperature = 0; // Kelvin
    simParams->drivingFreq = 32.25 * 1e9;
    simParams->oscillatingZeemanStrength = 1e-3;
    simParams->forceStopAtIteration = -1;
    simParams->gyroMagConst = GV.GetGyromagneticConstant();
    simParams->maxSimTime = 2e-9;
    simParams->satMag = -1;
    simParams->stepsize = 1e-15;

    // Data Output Parameters
    simStates->fixedOutputSites = {600, 1000, 1400, 1800, 2200, 2600, 3000, 3400};
    simParams->numberOfDataPoints = 1e2; //static_cast<int>(maxSimTime / recordingInterval);

    // Damping Factors
    simParams->gilbertDamping = 1e-4;
    simParams->gilbertABCInner = 1e-4;
    simParams->gilbertABCOuter = 1e0;

    // Spin chain and multi-layer Parameters
    simStates->discreteDrivenSites = {1};
    simParams->drivingRegionWidth = 57;
    simParams->numNeighbours = -1;
    simParams->numSpinsInABC = 300;
    simParams->numLayers = 1;

    // Shockwave Parameters (rarely used)
    simParams->iterStartShock = 0.0;
    simParams->iterEndShock = 0.0001;
    simParams->shockwaveGradientTime = 1;
    simParams->shockwaveInitialStrength = 0;  // Set equal to oscillatingZeemanStrength if NOT starting at time=0
    simParams->shockwaveMax = 3e-3;
    simParams->shockwaveScaling = 1;

    // Currently unused parameters
    simParams->recordingInterval = 2.5e-15;
    simParams->layerOfInterest = 1;

    // TODO Transferred from _generateRemainingParameters for now
    // Set flags from global variables
    simFlags->isFerromagnetic = GV.GetIsFerromagnetic();
    simParams->exchangeEnergyMin = GV.GetExchangeMinVal();
    simParams->exchangeEnergyMax = GV.GetExchangeMaxVal();
    simParams->dmiConstant = GV.GetDMIConstant() / 2; // Needed to match my Python and C++ code
    simParams->staticZeemanStrength = GV.GetStaticBiasField();
    simParams->anisotropyField = GV.GetAnisotropyField();
    simParams->numSpinsInChain = GV.GetNumSpins();
    simParams->exchangeStiffness = 5.3e-17;
    simParams->spinPolarisation = 0.5;
    simParams->spinTransferEfficiency = 0.4;
    simFlags->preferredDirection = 2;
    simParams->PERMEABILITY_IRON *= _BOHR_MAGNETON;  // Conversion to Am^2
    simParams->dipoleConstant = SimulationParameters::PERM_FREESPACE / (4.0 * M_PI);

    simParams->layerOfInterest -= 1;  // To correct for 0-indexing

    // Testing ONLY!
    simParams->numSpinsDRPeak = 57;
    simParams->numSpinsDRGradient =  (simParams->drivingRegionWidth - simParams->numSpinsDRPeak) / 2;

    // TODO. Turn into a struct
    simParams->dmiRegionLhs = -1;
    simParams->dmiRegionRhs = -1;
    simParams->dmiRegionOffset = 300;
    simParams->numSpinsDmiPeak = simParams->numSpinsDRPeak;
    simParams->numSpinsDmiWidth = -1; // This should be automatically calculated

    // TODO. Turn into a struct
    simParams->dampingRegionLhs = -1;
    simParams->dampingRegionRhs = -1;
    simParams->dampingRegionOffset = simParams->dmiRegionOffset;
    simParams->dampingGradientPeak = 1e-2;
    simParams->numSpinsDampingPeak = 657;
    simParams->numSpinsDampingWidth = -1; // This should be automatically calculated

}

void SolversInitialisation::_generateRemainingParameters() {

    // Computations based upon other inputs
    simParams->drivingAngFreq = 2 * M_PI * simParams->drivingFreq;

    simParams->iterationEnd = static_cast<int>(simParams->maxSimTime / simParams->stepsize);
    simParams->stepsizeHalf = simParams->stepsize / 2.0;

    simStates->layerSpinsInChain = {simParams->drivingRegionWidth, simParams->numSpinsInChain};
    simStates->layerSpinPairs.clear();
    simStates->layerTotalSpins.clear();
    for (int& spinsInChain: simStates->layerSpinsInChain) {
        simStates->layerSpinPairs.push_back(spinsInChain - 1);
        simStates->layerTotalSpins.push_back(spinsInChain + 2 * simParams->numSpinsInABC);
    }

    simStates->gilbertVectorMulti.resize(simParams->numLayers, {0});
    simParams->numberOfSpinPairs = simParams->numSpinsInChain - 1;

    //GV.SetNumSpins(simParams->numSpinsInChain + 2 * simParams->numSpinsInABC); // TODO turn this into a simParam
    simParams->systemTotalSpins = simParams->numSpinsInChain + 2 * simParams->numSpinsInABC;
}

void SolversInitialisation::_setMaterialParameters() {

    if (simFlags->isFerromagnetic)
        simParams->anisotropyField = 0;

    if (!simFlags->hasStaticZeeman)
        simParams->staticZeemanStrength = 0.0;
}

void SolversInitialisation::_guardClauses() {

    if (simFlags->shouldDriveCeaseEarly and simParams->iterEndShock <= 0) {
        std::cout << "Warning: [shouldDriveCeaseEarly: True] however [iterEndShock: " << simParams->iterEndShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (simFlags->hasShockwave and simParams->iterStartShock < 0) {
        std::cout << "Warning: [hasShockwave: True] however [iterStartShock: " << simParams->iterStartShock << " ! > 0.0]"
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

    if (simFlags->hasDemagIntense && simFlags->hasDemagFFT) {
        std::cout << "Warning: You cannot use both the intense and FFT demag solvers. Please choose one or the other.";
        exit(1);
    }

    if ((simFlags->hasDemagIntense && !GV.GetIsFerromagnetic()) || (simFlags->hasDemagFFT && !GV.GetIsFerromagnetic())) {
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

    if ( (simParams->numSpinsDRPeak + simParams->numSpinsDRGradient) > simParams->drivingRegionWidth) {
        std::cout << "Warning: Attempted DR gradient setup exceeds the drivingRegionWidth!";
        exit(1);
    }

}

void SolversInitialisation::testModifyingDouble(double newValue) {
    // Test modifying a double in the parent class from the child
    simParams->ambientTemperature = newValue;
    simParams->iterEndShock = 0.12345;
    std::cout << "I changed another value: " << simParams->iterEndShock << std::endl;
}