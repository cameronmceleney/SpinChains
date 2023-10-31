//
// Created by Cameron McEleney on 31/10/2023.
//

// C++ User Libraries (Parent)
#include "NMSuperClassTest.h"

// C++ User Libraries (Children)
#include "NMSubClasses/NMInitialisation.h"

// Define static variables
int NMSuperClassTest::iterationStart= 0;
double NMSuperClassTest::largestMNorm = 1e-50;
double NMSuperClassTest::PERMITTIVITY_IRON = 2.22; // cobalt = 1.72, iron = 2.22, nickel = 0.6;
double NMSuperClassTest::totalTime = 0;

// Test section
double              NMSuperClassTest::ambientTemperature = 0.0;
double              NMSuperClassTest::anisotropyField = 0.0;
double              NMSuperClassTest::drivingAngFreq = 0.0;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
double              NMSuperClassTest::drivingFreq = 0.0;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
int                 NMSuperClassTest::drivingRegionLhs = 0;                         // The position of the spin which is leftmost in the driving region.
int                 NMSuperClassTest::drivingRegionRhs = 0;                         // The position of the spin which is rightmost in the driving region.
int                 NMSuperClassTest::drivingRegionWidth;                       // Driving region width.
double              NMSuperClassTest::dynamicBiasField = 0.0;                         // Driving field amplitude [T] (caution: papers often give in [mT]).
int                 NMSuperClassTest::forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
double              NMSuperClassTest::dipoleConstant = 0.0;                           // Scaling factor which is constant across dipolar interaction calculations.
double              NMSuperClassTest::gilbertDamping = 0.0;                           // Gilbert damping factor for main chain.
double              NMSuperClassTest::gilbertABCInner = 0.0;                          // The lower Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the chain meets the ABC.
double              NMSuperClassTest::gilbertABCOuter = 0.0;                          // The upper Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the ABC meets the pinned sites.
double              NMSuperClassTest::gyroMagConst = 0.0;                             // Gyromagnetic ratio of an electron [GHz/T].
int                 NMSuperClassTest::iterationEnd = 0;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
double              NMSuperClassTest::iterStartShock = 0.0;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the maxSimTime.
double              NMSuperClassTest::iterEndShock = 0.0;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the maxSimTime.
double              NMSuperClassTest::maxSimTime = 0.0;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [mxInit + myInit + mzInit]  CANNOT sum to greater than 1.0
int                 NMSuperClassTest::numberNeighbours = 0;
int                 NMSuperClassTest::numberOfDataPoints = 0;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
int                 NMSuperClassTest::numberOfSpinPairs = 0;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
int                 NMSuperClassTest::numSpinsDamped = 0;                           // Number of spins in the damped regions (previously called _numGilbert).
int                 NMSuperClassTest::numSpinsInChain = 0;                          // The number of spin sites in the spin chain to be simulated.
int                 NMSuperClassTest::systemTotalSpins = 0;                         // The total number of spins in the system (chain plus ABCs).
double              NMSuperClassTest::satMag = 0.0;                                   // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)
double              NMSuperClassTest::shockwaveGradientTime = 0.0;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
double              NMSuperClassTest::shockwaveInitialStrength = 0.0;                 // Initial strength of the shockwave before shockwaveScaling occurs. (Default: = dynamicBiasField)
double              NMSuperClassTest::shockwaveMax = 0.0;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
double              NMSuperClassTest::shockwaveScaling = 0.0;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving
double              NMSuperClassTest::shockwaveStepsize = 0.0;                        // Size of incremental increase in shockwave amplitude.
double              NMSuperClassTest::stepsize = 0.0;                                 // Stepsize between values
double              NMSuperClassTest::stepsizeHalf = 0.0;                             // Separately defined to avoid repeated unnecessary calculations inside loops
std::string         NMSuperClassTest::stepsizeString;                           // Object to string conversation for stepsize
std::string         NMSuperClassTest::stopIterString;                           // Object to string conversion for _stopIterVal
int                 NMSuperClassTest::totalLayers = 0;


    // #####################################        Protected Flags        ####################################
bool NMSuperClassTest::centralDrive = false;                               // Drive from the centre of the chain if (true)
bool NMSuperClassTest::driveAllLayers = false;
bool NMSuperClassTest::dualDrive = false;
bool NMSuperClassTest::hasShockwave = false;
bool NMSuperClassTest::hasStaticDrive = false;
bool NMSuperClassTest::isShockwaveOn = false;
bool NMSuperClassTest::isShockwaveAtMax = false;
bool NMSuperClassTest::lhsDrive = false;
bool NMSuperClassTest::rhsDrive = false;
bool NMSuperClassTest::printAllData = false;
bool NMSuperClassTest::printFixedLines = false;
bool NMSuperClassTest::printFixedSites = false;
bool NMSuperClassTest::shouldDriveCease = false;
bool NMSuperClassTest::shouldTrackMValues = false;
bool NMSuperClassTest::useLLG = false;
bool NMSuperClassTest::useSLLG = false;
bool NMSuperClassTest::useDipolar = false;
bool NMSuperClassTest::useZeeman = false;
bool NMSuperClassTest::useDemagIntense = false;
bool NMSuperClassTest::useDemagFft = false;
bool NMSuperClassTest::useMultilayer = false;
bool NMSuperClassTest::debugFunc = false;

    // #####################################        Protected Data Structures        ####################################
std::vector<double> NMSuperClassTest::exchangeVec = {};
std::list <int> NMSuperClassTest::fixedOutputSites = {};
std::vector<int> NMSuperClassTest::layerSpinPairs = {};
std::vector<int>    NMSuperClassTest::layerSpinsInChain = {};
std::vector<int>    NMSuperClassTest::layerTotalSpins = {};
std::vector<std::vector<std::vector<double>>> NMSuperClassTest::m0Nest = {{{}}};
std::vector<std::vector<std::vector<double>>> NMSuperClassTest::m1Nest = {{{}}};
std::vector<std::vector<std::vector<double>>> NMSuperClassTest::m2Nest = {{{}}};



//
bool NMSuperClassTest::isFm = GV.GetIsFerromagnetic();
double NMSuperClassTest::exchangeEnergyMin = GV.GetExchangeMinVal();
double NMSuperClassTest::exchangeEnergyMax = GV.GetExchangeMaxVal();

std::vector<double> NMSuperClassTest::gilbertVector = {0};
std::vector<std::vector<double>> NMSuperClassTest::gilbertVectorMulti = {};
std::vector<double> NMSuperClassTest::largestMNormMulti = {1e-50, 1e-50};
std::vector<double> NMSuperClassTest::mx0 = {0};
std::vector<double> NMSuperClassTest::my0 = {0};
std::vector<double> NMSuperClassTest::mz0 = {0};

NMSuperClassTest::NMSuperClassTest() {
    NMInitialise = std::make_unique<NMInitialisation>();
    NMConfig = std::make_unique<NMConfiguration>();
}

// superMethod implementation
void NMSuperClassTest::execute() {
    NMInitialise->Initialise();
    NMConfig->Configure();
}

NMSuperClassTest::~NMSuperClassTest() {
    // Memory automatically deallocated by unique_ptr
}