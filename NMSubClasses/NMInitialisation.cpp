//
// Created by Cameron McEleney on 31/10/2023.
//

#include "NMInitialisation.h"

NMInitialisation::NMInitialisation(std::shared_ptr<SharedVariableHolder> data) : NMSuperClassTest(data) {}

void NMInitialisation::callInitialise() {
    // Virtual function access
    Initialise();
}

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
    sharedData -> shouldTrackMValues = true;

    // Model Type
    sharedData -> useLLG = true;
    sharedData -> useSLLG = false;

    // Interaction Flags
    sharedData -> hasShockwave = false;
    sharedData -> useDipolar = false;
    sharedData -> useZeeman = true;
    sharedData -> useDemagIntense = false;
    sharedData -> useDemagFft = false;

    // Material Flags
    sharedData -> useMultilayer = false;

    // Drive Flags
    sharedData -> centralDrive = false;
    sharedData -> driveAllLayers = false;
    sharedData -> dualDrive = false;
    sharedData -> lhsDrive = true;
    sharedData -> rhsDrive = false;
    sharedData -> hasStaticDrive = false;
    sharedData -> shouldDriveCease = false;

    // Output Flags
    sharedData -> printAllData = false;
    sharedData -> printFixedLines = true;
    sharedData -> printFixedSites = false;
}

void NMInitialisation::_setSimulationParameters() {

    // Main Parameters
    sharedData -> ambientTemperature = 273; // Kelvin
    std::cout << sharedData -> ambientTemperature << std::endl;

    sharedData -> drivingFreq = 43.5 * 1e9;
    sharedData -> dynamicBiasField = 3e-3;
    sharedData -> forceStopAtIteration = -1;
    sharedData -> gyroMagConst = GV.GetGyromagneticConstant();
    sharedData -> maxSimTime = 0.7e-9;
    sharedData -> satMag = 0.010032;
    sharedData -> stepsize = 1e-15;

    // Shockwave Parameters
    sharedData -> iterStartShock = 0.0;
    sharedData -> iterEndShock = 0.0001;
    sharedData -> shockwaveGradientTime = 1;
    sharedData -> shockwaveInitialStrength = 0;  // Set equal to dynamicBiasField if NOT starting at time=0
    sharedData -> shockwaveMax = 3e-3;
    sharedData -> shockwaveScaling = 1;

    // Data Output Parameters
    sharedData -> fixedOutputSites = {12158, 14529, 15320};
    // _recordingInterval = 1e-15;
    sharedData -> numberOfDataPoints = 100; //static_cast<int>(maxSimTime / recordingInterval);
    _layerOfInterest = 1;

    // Damping Factors
    sharedData -> gilbertABCInner = 1e-4;
    sharedData -> gilbertABCOuter = 1e0;

    // Spin chain and multi-layer Parameters
    sharedData -> drivingRegionWidth = 200;
    sharedData -> numberNeighbours = -1;
    sharedData -> numSpinsDamped = 0;
    sharedData -> totalLayers = 1;
}

void NMInitialisation::_generateRemainingParameters() {
    // Computations based upon other inputs
    sharedData -> drivingAngFreq = 2 * M_PI * sharedData -> drivingFreq;
    sharedData -> PERMITTIVITY_IRON *= _BOHR_MAGNETON;  // Conversion to Am^2
    sharedData -> dipoleConstant = SharedVariableHolder::PERM_FREESPACE / (4.0 * M_PI);

    sharedData -> iterationEnd = static_cast<int>(sharedData -> maxSimTime / sharedData -> stepsize);
    sharedData -> stepsizeHalf = sharedData -> stepsize / 2.0;

    sharedData -> numSpinsInChain = GV.GetNumSpins();
    sharedData -> numberOfSpinPairs = sharedData -> numSpinsInChain - 1;
    sharedData -> layerSpinsInChain = {sharedData -> drivingRegionWidth, sharedData -> numSpinsInChain};
    GV.SetNumSpins(sharedData -> numSpinsInChain + 2 * sharedData -> numSpinsDamped);
    sharedData -> systemTotalSpins = GV.GetNumSpins();

    sharedData -> layerSpinPairs.clear();
    sharedData -> layerTotalSpins.clear();
    for (int& spinsInChain: sharedData -> layerSpinsInChain) {
        sharedData -> layerSpinPairs.push_back(spinsInChain - 1);
        sharedData -> layerTotalSpins.push_back(spinsInChain + 2 * sharedData -> numSpinsDamped);
    }
    sharedData -> gilbertVectorMulti.resize(sharedData -> totalLayers, {0});

    _layerOfInterest -= 1;  // To correct for 0-indexing
}

void NMInitialisation::_setMaterialParameters() {

    if (sharedData -> isFm)
        sharedData -> anisotropyField = 0;
    else if (!sharedData -> isFm)
        sharedData -> anisotropyField = GV.GetAnisotropyField();

    if (!sharedData -> useZeeman)
        GV.SetStaticBiasField(0);
}

void NMInitialisation::_guardClauses() {

    if (sharedData -> shouldDriveCease and sharedData -> iterEndShock <= 0) {
        std::cout << "Warning: [shouldDriveCease: True] however [iterEndShock: " << sharedData -> iterEndShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (sharedData -> lhsDrive and sharedData -> rhsDrive) {
        std::cout << "Warning: [lhsDrive: True] and [rhsDrive: True] are both TRUE. Please choose one or the other."
                  << std::endl;
        exit(1);
    }

    if (sharedData -> hasShockwave and sharedData -> iterStartShock < 0) {
        std::cout << "Warning: [hasShockwave: True] however [iterStartShock: " << sharedData -> iterStartShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if ((sharedData -> printFixedSites and sharedData -> printFixedLines) or (sharedData -> printFixedSites and sharedData -> printAllData) or
        (sharedData -> printFixedLines and sharedData -> printAllData)) {
        std::cout << "Warning: Multiple output flags detected. [printFixedSites: " << sharedData -> printFixedSites
                  << "] | [printFixedLines: " << sharedData -> printFixedLines << "] | [printAllData: " << sharedData -> printAllData << "]"
                  << std::endl;
        exit(1);
    }

    if ((sharedData -> lhsDrive && sharedData -> centralDrive) || (sharedData -> lhsDrive && sharedData -> dualDrive) || (sharedData -> centralDrive && sharedData -> dualDrive)) {
        std::cout << "Warning: two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << sharedData -> lhsDrive << "\n_centralDrive: " << sharedData -> centralDrive << "\n_dualDrive: " << sharedData -> dualDrive
                  << "\n\nExiting...";
        exit(1);
    }

    if (sharedData -> printFixedSites and sharedData -> fixedOutputSites.empty()) {
        std::cout << "Warning: Request to print fixed sites, but no sites were given [fixedOutputSites: (";
        for (int & fixed_out_val : sharedData -> fixedOutputSites)
                std::cout << fixed_out_val << ", ";
        std::cout << ")].";
        exit(1);
    }

    if (sharedData -> numberOfDataPoints > sharedData -> iterationEnd) {
        std::cout << "Warning: You tried to print more data than was generated [numberOfDataPoints > iterationEnd]";
        exit(1);
    }

    if (sharedData -> useLLG and sharedData -> useSLLG) {
        std::cout << "Warning: You cannot use both the LLG and sLLG equations. Please choose one or the other.";
        exit(1);
    }

    if (sharedData -> useMultilayer and sharedData -> totalLayers < 2) {
        std::cout << "Warning: You cannot use the multilayer solver with less than 2 layers.";
        exit(1);
    }

    if (sharedData -> useDemagIntense && sharedData -> useDemagFft) {
        std::cout << "Warning: You cannot use both the intense and FFT demag solvers. Please choose one or the other.";
        exit(1);
    }

    if ((sharedData -> useDemagIntense && !GV.GetIsFerromagnetic()) || (sharedData -> useDemagFft && !GV.GetIsFerromagnetic())) {
        std::cout << "Warning: You cannot use the demag solvers with non-ferromagnetic materials.";
        exit(1);
    }

}
