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

    std::cout << "ambientTemp in init: " << systemData->ambientTemperature << std::endl;
    std::cout<<"Completed initialisation of init vals" << std::endl;
}
void NMInitialisation::_setSimulationFlags() {

    // Debugging Flags
    systemData->shouldTrackMValues = true;

    // Model Type
    systemData->useLLG = true;
    systemData->useSLLG = false;

    // Interaction Flags
    systemData->hasShockwave = false;
    systemData->useDipolar = false;
    systemData->useZeeman = true;
    systemData->useDemagIntense = false;
    systemData->useDemagFft = false;

    // Material Flags
    systemData->useMultilayer = false;

    // Drive Flags
    systemData->centralDrive = false;
    systemData->driveAllLayers = false;
    systemData->dualDrive = false;
    systemData->lhsDrive = true;
    systemData->rhsDrive = false;
    systemData->hasStaticDrive = false;
    systemData->shouldDriveCease = false;

    // Output Flags
    systemData->printAllData = false;
    systemData->printFixedLines = true;
    systemData->printFixedSites = false;
}

void NMInitialisation::_setSimulationParameters() {

    // Main Parameters
    systemData->ambientTemperature = 273; // Kelvin

    systemData->drivingFreq = 43.5 * 1e9;
    systemData->dynamicBiasField = 3e-3;
    systemData->forceStopAtIteration = -1;
    systemData->gyroMagConst = GV.GetGyromagneticConstant();
    systemData->maxSimTime = 0.7e-9;
    systemData->satMag = 0.010032;
    systemData->stepsize = 1e-15;

    // Shockwave Parameters
    systemData->iterStartShock = 0.0;
    systemData->iterEndShock = 0.0001;
    systemData->shockwaveGradientTime = 1;
    systemData->shockwaveInitialStrength = 0;  // Set equal to dynamicBiasField if NOT starting at time=0
    systemData->shockwaveMax = 3e-3;
    systemData->shockwaveScaling = 1;

    // Data Output Parameters
    systemData->fixedOutputSites = {12158, 14529, 15320};
    // _recordingInterval = 1e-15;
    systemData->numberOfDataPoints = 100; //static_cast<int>(maxSimTime / recordingInterval);
    _layerOfInterest = 1;

    // Damping Factors
    systemData->gilbertABCInner = 1e-4;
    systemData->gilbertABCOuter = 1e0;

    // Spin chain and multi-layer Parameters
    systemData->drivingRegionWidth = 200;
    systemData->numberNeighbours = -1;
    systemData->numSpinsDamped = 0;
    systemData->totalLayers = 1;
}

void NMInitialisation::_generateRemainingParameters() {
    // Computations based upon other inputs
    systemData->drivingAngFreq = 2 * M_PI * systemData->drivingFreq;
    systemData->PERMITTIVITY_IRON *= _BOHR_MAGNETON;  // Conversion to Am^2
    systemData->dipoleConstant = SystemDataContainer::PERM_FREESPACE / (4.0 * M_PI);

    systemData->iterationEnd = static_cast<int>(systemData->maxSimTime / systemData->stepsize);
    systemData->stepsizeHalf = systemData->stepsize / 2.0;

    systemData->numSpinsInChain = GV.GetNumSpins();
    systemData->numberOfSpinPairs = systemData->numSpinsInChain - 1;
    systemData->layerSpinsInChain = {systemData->drivingRegionWidth, systemData->numSpinsInChain};
    GV.SetNumSpins(systemData->numSpinsInChain + 2 * systemData->numSpinsDamped);
    systemData->systemTotalSpins = GV.GetNumSpins();

    systemData->layerSpinPairs.clear();
    systemData->layerTotalSpins.clear();
    for (int& spinsInChain: systemData->layerSpinsInChain) {
        systemData->layerSpinPairs.push_back(spinsInChain - 1);
        systemData->layerTotalSpins.push_back(spinsInChain + 2 * systemData->numSpinsDamped);
    }
    systemData->gilbertVectorMulti.resize(systemData->totalLayers, {0});

    _layerOfInterest -= 1;  // To correct for 0-indexing
}

void NMInitialisation::_setMaterialParameters() {

    if (systemData->isFm)
        systemData->anisotropyField = 0;
    else if (!systemData->isFm)
        systemData->anisotropyField = GV.GetAnisotropyField();

    if (!systemData->useZeeman)
        GV.SetStaticBiasField(0);
}

void NMInitialisation::_guardClauses() {

    if (systemData->shouldDriveCease and systemData->iterEndShock <= 0) {
        std::cout << "Warning: [shouldDriveCease: True] however [iterEndShock: " << systemData->iterEndShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (systemData->lhsDrive and systemData->rhsDrive) {
        std::cout << "Warning: [lhsDrive: True] and [rhsDrive: True] are both TRUE. Please choose one or the other."
                  << std::endl;
        exit(1);
    }

    if (systemData->hasShockwave and systemData->iterStartShock < 0) {
        std::cout << "Warning: [hasShockwave: True] however [iterStartShock: " << systemData->iterStartShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if ((systemData->printFixedSites and systemData->printFixedLines) or (systemData->printFixedSites and systemData->printAllData) or
        (systemData->printFixedLines and systemData->printAllData)) {
        std::cout << "Warning: Multiple output flags detected. [printFixedSites: " << systemData->printFixedSites
                  << "] | [printFixedLines: " << systemData->printFixedLines << "] | [printAllData: " << systemData->printAllData << "]"
                  << std::endl;
        exit(1);
    }

    if ((systemData->lhsDrive && systemData->centralDrive) || (systemData->lhsDrive && systemData->dualDrive) || (systemData->centralDrive && systemData->dualDrive)) {
        std::cout << "Warning: two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << systemData->lhsDrive << "\n_centralDrive: " << systemData->centralDrive << "\n_dualDrive: " << systemData->dualDrive
                  << "\n\nExiting...";
        exit(1);
    }

    if (systemData->printFixedSites and systemData->fixedOutputSites.empty()) {
        std::cout << "Warning: Request to print fixed sites, but no sites were given [fixedOutputSites: (";
        for (int & fixed_out_val : systemData->fixedOutputSites)
                std::cout << fixed_out_val << ", ";
        std::cout << ")].";
        exit(1);
    }

    if (systemData->numberOfDataPoints > systemData->iterationEnd) {
        std::cout << "Warning: You tried to print more data than was generated [numberOfDataPoints > iterationEnd]";
        exit(1);
    }

    if (systemData->useLLG and systemData->useSLLG) {
        std::cout << "Warning: You cannot use both the magDynamics and sLLG equations. Please choose one or the other.";
        exit(1);
    }

    if (systemData->useMultilayer and systemData->totalLayers < 2) {
        std::cout << "Warning: You cannot use the multilayer solver with less than 2 layers.";
        exit(1);
    }

    if (systemData->useDemagIntense && systemData->useDemagFft) {
        std::cout << "Warning: You cannot use both the intense and FFT demag solvers. Please choose one or the other.";
        exit(1);
    }

    if ((systemData->useDemagIntense && !GV.GetIsFerromagnetic()) || (systemData->useDemagFft && !GV.GetIsFerromagnetic())) {
        std::cout << "Warning: You cannot use the demag solvers with non-ferromagnetic materials.";
        exit(1);
    }

}

void NMInitialisation::testModifyingDouble(double newValue) {
    // Test modifying a double in the parent class from the child
    systemData->ambientTemperature = newValue;
    systemData->iterEndShock = 0.12345;
    std::cout << "I changed another value: " << systemData->iterEndShock << std::endl;
}