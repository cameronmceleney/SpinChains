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
}
void NMInitialisation::_setSimulationFlags() {

    // Debugging Flags
    shouldTrackMValues = true;

    // Model Type
    useLLG = true;
    useSLLG = false;

    // Interaction Flags
    hasShockwave = false;
    useDipolar = false;
    useZeeman = true;
    useDemagIntense = false;
    useDemagFft = false;

    // Material Flags
    useMultilayer = false;

    // Drive Flags
    centralDrive = false;
    driveAllLayers = false;
    dualDrive = false;
    lhsDrive = true;
    rhsDrive = false;
    hasStaticDrive = false;
    shouldDriveCease = false;

    // Output Flags
    printAllData = false;
    printFixedLines = true;
    printFixedSites = false;
}

void NMInitialisation::_setSimulationParameters() {

    // Main Parameters
    ambientTemperature = 273; // Kelvin
    drivingFreq = 43.5 * 1e9;
    dynamicBiasField = 3e-3;
    forceStopAtIteration = -1;
    gyroMagConst = GV.GetGyromagneticConstant();
    maxSimTime = 0.7e-9;
    satMag = 0.010032;
    stepsize = 1e-15;

    // Shockwave Parameters
    iterStartShock = 0.0;
    iterEndShock = 0.0001;
    shockwaveGradientTime = 1;
    shockwaveInitialStrength = 0;  // Set equal to dynamicBiasField if NOT starting at time=0
    shockwaveMax = 3e-3;
    shockwaveScaling = 1;

    // Data Output Parameters
    fixedOutputSites = {12158, 14529, 15320};
    // _recordingInterval = 1e-15;
    numberOfDataPoints = 100; //static_cast<int>(maxSimTime / recordingInterval);
    _layerOfInterest = 1;

    // Damping Factors
    gilbertABCInner = 1e-4;
    gilbertABCOuter = 1e0;

    // Spin chain and multi-layer Parameters
    drivingRegionWidth = 200;
    numberNeighbours = -1;
    numSpinsDamped = 0;
    totalLayers = 1;
}

void NMInitialisation::_generateRemainingParameters() {
    // Computations based upon other inputs
    drivingAngFreq = 2 * M_PI * drivingFreq;
    PERMITTIVITY_IRON *= _BOHR_MAGNETON;  // Conversion to Am^2
    dipoleConstant = PERM_FREESPACE / (4.0 * M_PI);

    iterationEnd = static_cast<int>(maxSimTime / stepsize);
    stepsizeHalf = stepsize / 2.0;

    numSpinsInChain = GV.GetNumSpins();
    numberOfSpinPairs = numSpinsInChain - 1;
    layerSpinsInChain = {drivingRegionWidth, numSpinsInChain};
    GV.SetNumSpins(numSpinsInChain + 2 * numSpinsDamped);
    systemTotalSpins = GV.GetNumSpins();

    layerSpinPairs.clear();
    layerTotalSpins.clear();
    for (int& spinsInChain: layerSpinsInChain) {
        layerSpinPairs.push_back(spinsInChain - 1);
        layerTotalSpins.push_back(spinsInChain + 2 * numSpinsDamped);
    }
    gilbertVectorMulti.resize(totalLayers, {0});

    _layerOfInterest -= 1;  // To correct for 0-indexing
}

void NMInitialisation::_setMaterialParameters() {

    if (isFm)
        anisotropyField = 0;
    else if (!isFm)
        anisotropyField = GV.GetAnisotropyField();

    if (!useZeeman)
        GV.SetStaticBiasField(0);
}

void NMInitialisation::_guardClauses() {

    if (shouldDriveCease and iterEndShock <= 0) {
        std::cout << "Warning: [shouldDriveCease: True] however [iterEndShock: " << iterEndShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (lhsDrive and rhsDrive) {
        std::cout << "Warning: [lhsDrive: True] and [rhsDrive: True] are both TRUE. Please choose one or the other."
                  << std::endl;
        exit(1);
    }

    if (hasShockwave and iterStartShock < 0) {
        std::cout << "Warning: [hasShockwave: True] however [iterStartShock: " << iterStartShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if ((printFixedSites and printFixedLines) or (printFixedSites and printAllData) or
        (printFixedLines and printAllData)) {
        std::cout << "Warning: Multiple output flags detected. [printFixedSites: " << printFixedSites
                  << "] | [printFixedLines: " << printFixedLines << "] | [printAllData: " << printAllData << "]"
                  << std::endl;
        exit(1);
    }

    if ((lhsDrive && centralDrive) || (lhsDrive && dualDrive) || (centralDrive && dualDrive)) {
        std::cout << "Warning: two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << lhsDrive << "\n_centralDrive: " << centralDrive << "\n_dualDrive: " << dualDrive
                  << "\n\nExiting...";
        exit(1);
    }

    if (printFixedSites and fixedOutputSites.empty()) {
        std::cout << "Warning: Request to print fixed sites, but no sites were given [fixedOutputSites: (";
        for (int & fixed_out_val : fixedOutputSites)
                std::cout << fixed_out_val << ", ";
        std::cout << ")].";
        exit(1);
    }

    if (numberOfDataPoints > iterationEnd) {
        std::cout << "Warning: You tried to print more data than was generated [numberOfDataPoints > iterationEnd]";
        exit(1);
    }

    if (useLLG and useSLLG) {
        std::cout << "Warning: You cannot use both the LLG and sLLG equations. Please choose one or the other.";
        exit(1);
    }

    if (useMultilayer and totalLayers < 2) {
        std::cout << "Warning: You cannot use the multilayer solver with less than 2 layers.";
        exit(1);
    }

    if (useDemagIntense && useDemagFft) {
        std::cout << "Warning: You cannot use both the intense and FFT demag solvers. Please choose one or the other.";
        exit(1);
    }

    if ((useDemagIntense && !GV.GetIsFerromagnetic()) || (useDemagFft && !GV.GetIsFerromagnetic())) {
        std::cout << "Warning: You cannot use the demag solvers with non-ferromagnetic materials.";
        exit(1);
    }

}
