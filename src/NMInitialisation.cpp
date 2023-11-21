//
// Created by Cameron McEleney on 31/10/2023.
//

#include <utility>

#include "../include/NMInitialisation.h"

NMInitialisation::NMInitialisation(std::shared_ptr<SimulationParameters> paramsData,
                                   std::shared_ptr<SimulationStates> sharedSimStates,
                                   std::shared_ptr<SimulationFlags> sharedSimFlags)  :

   SolversSuperClass(std::move(paramsData), std::move(sharedSimStates), std::move(sharedSimFlags)) {}

void NMInitialisation::Initialise() {
    // Constructor implementation
    _setSimulationFlags();
    _setSimulationParameters();
    _generateRemainingParameters();
    _setMaterialParameters();
    _guardClauses();
}
void NMInitialisation::_setSimulationFlags() {

    // Debugging Flags
    simFlags->shouldTrackMValues = true;

    // Model Type
    simFlags->useLLG = true;
    simFlags->useSLLG = false;

    // Interaction Flags
    simFlags->hasShockwave = false;
    simFlags->useDipolar = false;
    simFlags->useZeeman = true;
    simFlags->useDemagIntense = false;
    simFlags->useDemagFft = false;

    // Material Flags
    simFlags->useMultilayer = false;

    // Drive Flags
    simFlags->centralDrive = false;
    simFlags->driveAllLayers = false;
    simFlags->dualDrive = false;
    simFlags->lhsDrive = false;
    simFlags->rhsDrive = true;
    simFlags->hasStaticDrive = false;
    simFlags->shouldDriveCease = false;

    // Output Flags
    simFlags->printAllData = false;
    simFlags->printFixedLines = true;
    simFlags->printFixedSites = false;
}

void NMInitialisation::_setSimulationParameters() {


    simFlags->isFm = GV.GetIsFerromagnetic();
    simParams->exchangeEnergyMin = GV.GetExchangeMinVal();
    simParams->exchangeEnergyMax = GV.GetExchangeMaxVal();

    // Main Parameters
    simParams->ambientTemperature = 273; // Kelvin

    simParams->drivingFreq = 43.5 * 1e9;
    simParams->dynamicBiasField = 3e-3;
    simParams->forceStopAtIteration = -1;
    simParams->gyroMagConst = GV.GetGyromagneticConstant();
    simParams->maxSimTime = 0.7e-9;
    simParams->satMag = 0.010032;
    simParams->stepsize = 1e-15;

    // Shockwave Parameters
    simParams->iterStartShock = 0.0;
    simParams->iterEndShock = 0.0001;
    simParams->shockwaveGradientTime = 1;
    simParams->shockwaveInitialStrength = 0;  // Set equal to dynamicBiasField if NOT starting at time=0
    simParams->shockwaveMax = 3e-3;
    simParams->shockwaveScaling = 1;

    // Data Output Parameters
    simParams->fixedOutputSites = {12158, 14529, 15320};
    // _recordingInterval = 1e-15;
    simParams->numberOfDataPoints = 100; //static_cast<int>(maxSimTime / recordingInterval);
    _layerOfInterest = 1;

    // Damping Factors
    simParams->gilbertABCInner = 1e-4;
    simParams->gilbertABCOuter = 1e0;

    // Spin chain and multi-layer Parameters
    simParams->drivingRegionWidth = 200;
    simParams->numberNeighbours = -1;
    simParams->numSpinsDamped = 0;
    simParams->totalLayers = 1;
}

void NMInitialisation::_generateRemainingParameters() {
    // Computations based upon other inputs
    simParams->drivingAngFreq = 2 * M_PI * simParams->drivingFreq;
    simParams->PERMITTIVITY_IRON *= _BOHR_MAGNETON;  // Conversion to Am^2
    simParams->dipoleConstant = SimulationParameters::PERM_FREESPACE / (4.0 * M_PI);

    simParams->iterationEnd = static_cast<int>(simParams->maxSimTime / simParams->stepsize);
    simParams->stepsizeHalf = simParams->stepsize / 2.0;

    simParams->numSpinsInChain = GV.GetNumSpins();
    simParams->numberOfSpinPairs = simParams->numSpinsInChain - 1;
    simStates->layerSpinsInChain = {simParams->drivingRegionWidth, simParams->numSpinsInChain};
    GV.SetNumSpins(simParams->numSpinsInChain + 2 * simParams->numSpinsDamped);
    simParams->systemTotalSpins = GV.GetNumSpins();

    simStates->layerSpinPairs.clear();
    simStates->layerTotalSpins.clear();
    for (int& spinsInChain: simStates->layerSpinsInChain) {
        simStates->layerSpinPairs.push_back(spinsInChain - 1);
        simStates->layerTotalSpins.push_back(spinsInChain + 2 * simParams->numSpinsDamped);
    }
    simStates->gilbertVectorMulti.resize(simParams->totalLayers, {0});

    _layerOfInterest -= 1;  // To correct for 0-indexing
}

void NMInitialisation::_setMaterialParameters() {

    if (simFlags->isFm)
        simParams->anisotropyField = 0;
    else if (!simFlags->isFm)
        simParams->anisotropyField = GV.GetAnisotropyField();

    if (!simFlags->useZeeman)
        GV.SetStaticBiasField(0);
}

void NMInitialisation::_guardClauses() {

    if (simFlags->shouldDriveCease and simParams->iterEndShock <= 0) {
        std::cout << "Warning: [shouldDriveCease: True] however [iterEndShock: " << simParams->iterEndShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (simFlags->lhsDrive and simFlags->rhsDrive) {
        std::cout << "Warning: [lhsDrive: True] and [rhsDrive: True] are both TRUE. Please choose one or the other."
                  << std::endl;
        exit(1);
    }

    if (simFlags->hasShockwave and simParams->iterStartShock < 0) {
        std::cout << "Warning: [hasShockwave: True] however [iterStartShock: " << simParams->iterStartShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if ((simFlags->printFixedSites and simFlags->printFixedLines) or (simFlags->printFixedSites and simFlags->printAllData) or
        (simFlags->printFixedLines and simFlags->printAllData)) {
        std::cout << "Warning: Multiple output flags detected. [printFixedSites: " << simFlags->printFixedSites
                  << "] | [printFixedLines: " << simFlags->printFixedLines << "] | [printAllData: " << simFlags->printAllData << "]"
                  << std::endl;
        exit(1);
    }

    if ((simFlags->lhsDrive && simFlags->centralDrive) || (simFlags->lhsDrive && simFlags->dualDrive) || (simFlags->centralDrive && simFlags->dualDrive)) {
        std::cout << "Warning: two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << simFlags->lhsDrive << "\n_centralDrive: " << simFlags->centralDrive << "\n_dualDrive: " << simFlags->dualDrive
                  << "\n\nExiting...";
        exit(1);
    }

    if (simFlags->printFixedSites and simParams->fixedOutputSites.empty()) {
        std::cout << "Warning: Request to print fixed sites, but no sites were given [fixedOutputSites: (";
        for (int & fixed_out_val : simParams->fixedOutputSites)
                std::cout << fixed_out_val << ", ";
        std::cout << ")].";
        exit(1);
    }

    if (simParams->numberOfDataPoints > simParams->iterationEnd) {
        std::cout << "Warning: You tried to print more data than was generated [numberOfDataPoints > iterationEnd]";
        exit(1);
    }

    if (simFlags->useLLG and simFlags->useSLLG) {
        std::cout << "Warning: You cannot use both the magDynamics and sLLG equations. Please choose one or the other.";
        exit(1);
    }

    if (simFlags->useMultilayer and simParams->totalLayers < 2) {
        std::cout << "Warning: You cannot use the multilayer solver with less than 2 layers.";
        exit(1);
    }

    if (simFlags->useDemagIntense && simFlags->useDemagFft) {
        std::cout << "Warning: You cannot use both the intense and FFT demag solvers. Please choose one or the other.";
        exit(1);
    }

    if ((simFlags->useDemagIntense && !GV.GetIsFerromagnetic()) || (simFlags->useDemagFft && !GV.GetIsFerromagnetic())) {
        std::cout << "Warning: You cannot use the demag solvers with non-ferromagnetic materials.";
        exit(1);
    }

}

void NMInitialisation::testModifyingDouble(double newValue) {
    // Test modifying a double in the parent class from the child
    simParams->ambientTemperature = newValue;
    simParams->iterEndShock = 0.12345;
    std::cout << "I changed another value: " << simParams->iterEndShock << std::endl;
}