//
// Created by Cameron McEleney on 31/10/2023.
//

#include "NMInitialisation.h"

void NMInitialisation::Initialise() {
    // Constructor implementation
    _setSimulationFlags();
    _setSimulationParameters();
    _generateRemainingParameters();
    _setMaterialParameters();
    _guardClauses();

    std::cout << "ambientTemp in init: " << simState->ambientTemperature << std::endl;
    std::cout<<"Completed initialisation of init vals" << std::endl;
}
void NMInitialisation::_setSimulationFlags() {

    // Debugging Flags
    simState->shouldTrackMValues = true;

    // Model Type
    simState->useLLG = true;
    simState->useSLLG = false;

    // Interaction Flags
    simState->hasShockwave = false;
    simState->useDipolar = false;
    simState->useZeeman = true;
    simState->useDemagIntense = false;
    simState->useDemagFft = false;

    // Material Flags
    simState->useMultilayer = false;

    // Drive Flags
    simState->centralDrive = false;
    simState->driveAllLayers = false;
    simState->dualDrive = false;
    simState->lhsDrive = false;
    simState->rhsDrive = true;
    simState->hasStaticDrive = false;
    simState->shouldDriveCease = false;

    // Output Flags
    simState->printAllData = false;
    simState->printFixedLines = true;
    simState->printFixedSites = false;
}

void NMInitialisation::_setSimulationParameters() {


    simState->isFm = GV.GetIsFerromagnetic();
    simState->exchangeEnergyMin = GV.GetExchangeMinVal();
    simState->exchangeEnergyMax = GV.GetExchangeMaxVal();

    // Main Parameters
    simState->ambientTemperature = 273; // Kelvin

    simState->drivingFreq = 43.5 * 1e9;
    simState->dynamicBiasField = 3e-3;
    simState->forceStopAtIteration = -1;
    simState->gyroMagConst = GV.GetGyromagneticConstant();
    simState->maxSimTime = 0.7e-9;
    simState->satMag = 0.010032;
    simState->stepsize = 1e-15;

    // Shockwave Parameters
    simState->iterStartShock = 0.0;
    simState->iterEndShock = 0.0001;
    simState->shockwaveGradientTime = 1;
    simState->shockwaveInitialStrength = 0;  // Set equal to dynamicBiasField if NOT starting at time=0
    simState->shockwaveMax = 3e-3;
    simState->shockwaveScaling = 1;

    // Data Output Parameters
    simState->fixedOutputSites = {12158, 14529, 15320};
    // _recordingInterval = 1e-15;
    simState->numberOfDataPoints = 100; //static_cast<int>(maxSimTime / recordingInterval);
    _layerOfInterest = 1;

    // Damping Factors
    simState->gilbertABCInner = 1e-4;
    simState->gilbertABCOuter = 1e0;

    // Spin chain and multi-layer Parameters
    simState->drivingRegionWidth = 200;
    simState->numberNeighbours = -1;
    simState->numSpinsDamped = 0;
    simState->totalLayers = 1;
}

void NMInitialisation::_generateRemainingParameters() {
    // Computations based upon other inputs
    simState->drivingAngFreq = 2 * M_PI * simState->drivingFreq;
    simState->PERMITTIVITY_IRON *= _BOHR_MAGNETON;  // Conversion to Am^2
    simState->dipoleConstant = SystemDataContainer::PERM_FREESPACE / (4.0 * M_PI);

    simState->iterationEnd = static_cast<int>(simState->maxSimTime / simState->stepsize);
    simState->stepsizeHalf = simState->stepsize / 2.0;

    simState->numSpinsInChain = GV.GetNumSpins();
    simState->numberOfSpinPairs = simState->numSpinsInChain - 1;
    simState->layerSpinsInChain = {simState->drivingRegionWidth, simState->numSpinsInChain};
    GV.SetNumSpins(simState->numSpinsInChain + 2 * simState->numSpinsDamped);
    simState->systemTotalSpins = GV.GetNumSpins();

    simState->layerSpinPairs.clear();
    simState->layerTotalSpins.clear();
    for (int& spinsInChain: simState->layerSpinsInChain) {
        simState->layerSpinPairs.push_back(spinsInChain - 1);
        simState->layerTotalSpins.push_back(spinsInChain + 2 * simState->numSpinsDamped);
    }
    simState->gilbertVectorMulti.resize(simState->totalLayers, {0});

    _layerOfInterest -= 1;  // To correct for 0-indexing
}

void NMInitialisation::_setMaterialParameters() {

    if (simState->isFm)
        simState->anisotropyField = 0;
    else if (!simState->isFm)
        simState->anisotropyField = GV.GetAnisotropyField();

    if (!simState->useZeeman)
        GV.SetStaticBiasField(0);
}

void NMInitialisation::_guardClauses() {

    if (simState->shouldDriveCease and simState->iterEndShock <= 0) {
        std::cout << "Warning: [shouldDriveCease: True] however [iterEndShock: " << simState->iterEndShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (simState->lhsDrive and simState->rhsDrive) {
        std::cout << "Warning: [lhsDrive: True] and [rhsDrive: True] are both TRUE. Please choose one or the other."
                  << std::endl;
        exit(1);
    }

    if (simState->hasShockwave and simState->iterStartShock < 0) {
        std::cout << "Warning: [hasShockwave: True] however [iterStartShock: " << simState->iterStartShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if ((simState->printFixedSites and simState->printFixedLines) or (simState->printFixedSites and simState->printAllData) or
        (simState->printFixedLines and simState->printAllData)) {
        std::cout << "Warning: Multiple output flags detected. [printFixedSites: " << simState->printFixedSites
                  << "] | [printFixedLines: " << simState->printFixedLines << "] | [printAllData: " << simState->printAllData << "]"
                  << std::endl;
        exit(1);
    }

    if ((simState->lhsDrive && simState->centralDrive) || (simState->lhsDrive && simState->dualDrive) || (simState->centralDrive && simState->dualDrive)) {
        std::cout << "Warning: two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << simState->lhsDrive << "\n_centralDrive: " << simState->centralDrive << "\n_dualDrive: " << simState->dualDrive
                  << "\n\nExiting...";
        exit(1);
    }

    if (simState->printFixedSites and simState->fixedOutputSites.empty()) {
        std::cout << "Warning: Request to print fixed sites, but no sites were given [fixedOutputSites: (";
        for (int & fixed_out_val : simState->fixedOutputSites)
                std::cout << fixed_out_val << ", ";
        std::cout << ")].";
        exit(1);
    }

    if (simState->numberOfDataPoints > simState->iterationEnd) {
        std::cout << "Warning: You tried to print more data than was generated [numberOfDataPoints > iterationEnd]";
        exit(1);
    }

    if (simState->useLLG and simState->useSLLG) {
        std::cout << "Warning: You cannot use both the magDynamics and sLLG equations. Please choose one or the other.";
        exit(1);
    }

    if (simState->useMultilayer and simState->totalLayers < 2) {
        std::cout << "Warning: You cannot use the multilayer solver with less than 2 layers.";
        exit(1);
    }

    if (simState->useDemagIntense && simState->useDemagFft) {
        std::cout << "Warning: You cannot use both the intense and FFT demag solvers. Please choose one or the other.";
        exit(1);
    }

    if ((simState->useDemagIntense && !GV.GetIsFerromagnetic()) || (simState->useDemagFft && !GV.GetIsFerromagnetic())) {
        std::cout << "Warning: You cannot use the demag solvers with non-ferromagnetic materials.";
        exit(1);
    }

}

void NMInitialisation::testModifyingDouble(double newValue) {
    // Test modifying a double in the parent class from the child
    simState->ambientTemperature = newValue;
    simState->iterEndShock = 0.12345;
    std::cout << "I changed another value: " << simState->iterEndShock << std::endl;
}