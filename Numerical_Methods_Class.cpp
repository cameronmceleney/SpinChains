#include "Numerical_Methods_Class.h"

void Numerical_Methods_Class::NMSetup() {
    // Core Flags
    _hasShockwave = false;
    _isFM = GV.GetIsFerromagnetic();
    _shouldTrackMValues = true;
    _useLLG = true;
    _useDipolar = false;
    _useBilayer = false;

    // Drive Flags
    _centralDrive = false;
    _dualDrive = false;
    _lhsDrive = true;
    _hasStaticDrive = false;
    _shouldDriveCease = false;

    // Core Parameters
    _recordingInterval = .7e-11;
    _drivingFreq = 42.5 * 1e9;
    _dynamicBiasField = 3e-3;
    _forceStopAtIteration = -1;
    _gyroMagConst = GV.GetGyromagneticConstant();
    _maxSimTime = 0.7e-9;
    _stepsize = 1e-15;

    // Shockwave Parameters
    _iterStartShock = 0.0;
    _iterEndShock = 0.0001;
    _shockwaveGradientTime = 1;
    _shockwaveInitialStrength = 0;  // Set equal to _dynamicBiasField if NOT starting at time=0
    _shockwaveMax = 3e-3;
    _shockwaveScaling = 1;

    // Data Output Parameters
    _fixed_output_sites = {12158, 14529, 15320};
    _numberOfDataPoints = static_cast<int>(_maxSimTime / _recordingInterval);

    _printAllData = false;
    _printFixedLines = true;
    _printFixedSites = false;

    // Damping Factors
    _gilbertConst  = 1e-4;
    _gilbertLower = _gilbertConst;
    _gilbertUpper = 1e0;

    // SpinChain Length Parameters
    _drivingRegionWidth = 200;
    _numSpinsDamped = 0;

    // Computations based upon other inputs
    _drivingAngFreq = 2 * M_PI * _drivingFreq;
    _iterationEnd = static_cast<int>(_maxSimTime / _stepsize);
    _numSpinsInChain = GV.GetNumSpins();
    _numberOfSpinPairs = _numSpinsInChain - 1;
    GV.SetNumSpins(_numSpinsInChain + 2 * _numSpinsDamped);
    _stepsizeHalf = _stepsize / 2.0;

    if (_isFM)
        _anisotropyField = 0;
    else if (!_isFM)
        _anisotropyField = GV.GetAnisotropyField();

    // ###################### Core Method Invocations ######################
    // Order is intentional, and must be maintained!
    FinalChecks();
    SetShockwaveConditions();
    SetDampingRegion();
    SetDrivingRegion();
    SetExchangeVector();
    SetInitialMagneticMoments();
}
void Numerical_Methods_Class::SetDampingRegion() {
    // Generate the damping regions that are appended to either end of the spin chain.

    LinspaceClass DampingRegionLeft;
    LinspaceClass DampingRegionRight;

    if (_numSpinsDamped < 0) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(0);
    }

    std::vector<double> gilbertChain(_numSpinsInChain, _gilbertConst);

    DampingRegionLeft.set_values(_gilbertUpper, _gilbertLower, _numSpinsDamped, true, false);
    DampingRegionRight.set_values(_gilbertLower, _gilbertUpper, _numSpinsDamped, true, false);
    std::vector<double> tempGilbertLHS = DampingRegionLeft.generate_array();
    std::vector<double> tempGilbertRHS = DampingRegionRight.generate_array();

    // Combine all damped regions to form vector which describes the entire spinchain.
    _gilbertVector.insert(_gilbertVector.end(), tempGilbertLHS.begin(), tempGilbertLHS.end());
    _gilbertVector.insert(_gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());
    _gilbertVector.insert(_gilbertVector.end(), tempGilbertRHS.begin(), tempGilbertRHS.end());
    _gilbertVector.push_back(0);
}
void Numerical_Methods_Class::SetDrivingRegion() {
    /**
     * Set up driving regions for the system. The LHS option is solely for drives from the left of the system. The RHS options contains the
     * drive from the right, as well as an option to drive from the centre.
     */

    if (_centralDrive) {
        _drivingRegionLHS = (_numSpinsInChain/2) +_numSpinsDamped - (_drivingRegionWidth / 2);
        _drivingRegionRHS = (_numSpinsInChain/2) +_numSpinsDamped + (_drivingRegionWidth / 2);
        return;
    }

    if (_dualDrive) {
        return;
    }

    if (_lhsDrive) {
        // The +1/-1 offset excludes the zeroth spin while retaining the correct driving width
        _drivingRegionLHS = _numSpinsDamped + 1;
        _drivingRegionRHS = _drivingRegionLHS + _drivingRegionWidth - 1;
        return;
    }

    if (!_lhsDrive) {
        // The +1 is to correct the offset of adding a zeroth spin
        _drivingRegionRHS = GV.GetNumSpins() - _numSpinsDamped - 1;
        _drivingRegionLHS = _drivingRegionRHS - _drivingRegionWidth + 1;
        return;
    }

}
void Numerical_Methods_Class::SetExchangeVector() {
    /*
     * Create the arrays which house the exchange integral values. There are options to have a non-uniform exchange coded in, as well as the option to
     * induce a 'kick' into the system by initialising certain spins to have differing parameters to their neighbours.
     */
    LinspaceClass SpinChainExchange;

    if (_numSpinsDamped > 0) {
        SpinChainExchange.set_values(GV.GetExchangeMinVal(), GV.GetExchangeMaxVal(), _numberOfSpinPairs, true, false);
        _exchangeVec = SpinChainExchange.generate_array();

        std::vector<double> dampingRegionLeftExchange(_numSpinsDamped, GV.GetExchangeMinVal()), dampingRegionRightExchange(_numSpinsDamped, GV.GetExchangeMaxVal());
        dampingRegionLeftExchange.insert(dampingRegionLeftExchange.begin(), 0);
        dampingRegionRightExchange.push_back(0);

        _exchangeVec.insert(_exchangeVec.begin(), dampingRegionLeftExchange.begin(), dampingRegionLeftExchange.end());
        _exchangeVec.insert(_exchangeVec.end(), dampingRegionRightExchange.begin(), dampingRegionRightExchange.end());
    } else {
        // The linearly spaced vector is saved as the class member '_exchangeVec' simply to increase code readability
        SpinChainExchange.set_values(GV.GetExchangeMinVal(), GV.GetExchangeMaxVal(), _numberOfSpinPairs, true, true);
        _exchangeVec = SpinChainExchange.generate_array();
    }
}
void Numerical_Methods_Class::FinalChecks() {

    if (_shouldDriveCease and _iterEndShock <= 0) {
        std::cout << "Warning: [_shouldDriveCease: True] however [_iterEndShock: " << _iterEndShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if (_hasShockwave and _iterStartShock < 0) {
        std::cout << "Warning: [_hasShockwave: True] however [_iterStartShock: " << _iterStartShock << " ! > 0.0]"
                  << std::endl;
        exit(1);
    }

    if ((_printFixedSites and _printFixedLines) or (_printFixedSites and _printAllData) or
        (_printFixedLines and _printAllData)) {
        std::cout << "Warning: Multiple output flags detected. [_printFixedSites: " << _printFixedSites
                  << "] | [_printFixedLines: " << _printFixedLines << "] | [_printAllData: " << _printAllData << "]"
                  << std::endl;
        exit(1);
    }

    if ((_lhsDrive && _centralDrive) || (_lhsDrive && _dualDrive) || (_centralDrive && _dualDrive)) {
        std::cout << "Warning: two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << _lhsDrive << "\n_centralDrive: " << _centralDrive << "\n_dualDrive: " << _dualDrive
                  << "\n\nExiting...";
        exit(1);
    }

    if (_printFixedSites and _fixed_output_sites.empty()) {
        std::cout << "Warning: Request to print fixed sites, but no sites were given [_fixed_output_sites: (";
        for (int & fixed_out_val : _fixed_output_sites)
                std::cout << fixed_out_val << ", ";
        std::cout << ")].";
        exit(1);
    }

    if (_numberOfDataPoints > _iterationEnd) {
        std::cout << "Warning: You tried to print more data than was generated [_numberOfDataPoints > _iterationEnd]";
        exit(1);
    }
}
void Numerical_Methods_Class::SetInitialMagneticMoments() {

    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mxInitCond(GV.GetNumSpins(), _mxInit), myInitCond(GV.GetNumSpins(), _myInit), mzInitCond(GV.GetNumSpins(), _mzInit);
    // mxInitCond[0] = _mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < GV.GetNumSpins(); i++) {
        mxInitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    if (!_isFM) {
        for (int i = 0; i < GV.GetNumSpins(); i++) {
            if (i % 2 == 1)
                mzInitCond[i] *= -1.0;
        }
    }

    // Appends initial conditions to the vectors
    _mx0.insert(_mx0.end(), mxInitCond.begin(), mxInitCond.end());
    _my0.insert(_my0.end(), myInitCond.begin(), myInitCond.end());
    _mz0.insert(_mz0.end(), mzInitCond.begin(), mzInitCond.end());

    // This zero is the (N+1)th spin on the RHS of the chain
    _mx0.push_back(0);
    _my0.push_back(0);
    _mz0.push_back(0);
}
void Numerical_Methods_Class::SetShockwaveConditions() {

    if (_hasShockwave) {
        _shockwaveStepsize = (_shockwaveMax - _shockwaveInitialStrength) / _shockwaveGradientTime;
    } else {
        // Ensures, on the output file, all parameter read as zero; reduces confusion when no shockwave is applied.
        _iterStartShock = 0;
        _shockwaveScaling = 0;
        _shockwaveGradientTime = 0;
        _shockwaveInitialStrength = 0;
        _shockwaveMax = _shockwaveInitialStrength * _shockwaveScaling;
        _shockwaveStepsize = (_shockwaveMax - _shockwaveInitialStrength) / _shockwaveGradientTime;
    }
}

void Numerical_Methods_Class::SetInitialMagneticMomentsMultilayer(std::vector<std::vector<std::vector<double>>>& nestedNestedVector,
                                                                  int layer, double mxInit, double myInit, double mzInit) {
    layer -= 1;
    // mxInitCond[0] = _mxInit; // Only perturb initial spin

    /*
    for (int i = 0; i < GV.GetNumSpins(); i++) {
        mxInitCond[i] = 0.003162277;
        // myInitCond[i] = 0.0;
        mzInitCond[i] = 0.999994999;
    }
    */

    for (int i = 0; i < GV.GetNumSpins(); i++) {
        nestedNestedVector[layer].push_back({mxInit, myInit, mzInit});
    }

    // This zero is the (N+1)th spin on the RHS of the chain
    nestedNestedVector[layer].push_back({0.0, 0.0, 0.0});

}
std::vector<std::vector<std::vector<double>>> Numerical_Methods_Class::initializeNestedNestedVector(int numLayers, bool includeEnd) {
    std::vector<std::vector<std::vector<double>>> innerNestedVector;
    for (int j = 0; j < numLayers; j++) {
        std::vector<std::vector<double>> innerVector;
        std::vector<double> innermostVector = {0.0, 0.0, 0.0};
        innerVector.push_back(innermostVector);
        innerNestedVector.push_back(innerVector);
    }
    return innerNestedVector;
}


std::vector<double> Numerical_Methods_Class::DipoleDipoleCoupling(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                                  std::vector<double> mzTerms, std::vector<int> sitePositions) {
    std::vector<double> totalDipoleTerms = {0.0, 0.0, 0.0};
    //double mu1x = mxTerms[1];
    //double mu1y = myTerms[1];
    //double mu1z = mzTerms[1];

    double exchangeStiffness = 5.3e-17;

    for (int i = 0; i < mxTerms.size(); i++) {

        if (i == 1) {
            continue;
        }

        double latticeConstant = sqrt(exchangeStiffness / _exchangeVec[i]);

        std::vector<double> positionVector = {(sitePositions[i] - sitePositions[1]) * latticeConstant, 0, 0};

        //std::vector<double> positionVector = {  (mxTerms[i] -  mu1x) * latticeConstant,
        //                                        (myTerms[i] -  mu1y) * latticeConstant,
        //                                        (mzTerms[i] -  mu1z) * latticeConstant};

        double positionVector_norm = std::sqrt(positionVector[0] * positionVector[0] + positionVector[1] * positionVector[1] + positionVector[2] * positionVector[2]);
        double positionVector_cubed = std::pow(positionVector_norm, 3);
        double positionVector_fifth = std::pow(positionVector_norm, 5);

        double mu1DotPosition = mxTerms[1] * positionVector[0] + myTerms[1] * positionVector[1] + mzTerms[1] * positionVector[2];
        double mu2DotPosition = mxTerms[i] * positionVector[0] + myTerms[i] * positionVector[1] + mzTerms[i] * positionVector[2];
        double mu1_dot_mu2 = mxTerms[1] * mxTerms[i];

        double gConstant = _permFreeSpace / (4.0 * M_PI);
        std::vector<double> DipoleValues = {gConstant* ((3 * positionVector[0] * mu2DotPosition)/positionVector_fifth - mxTerms[i]/positionVector_cubed),
                                            gConstant* ((3 * positionVector[1] * mu2DotPosition)/positionVector_fifth - myTerms[i]/positionVector_cubed),
                                            gConstant* ((3 * positionVector[2] * mu2DotPosition)/positionVector_fifth - mzTerms[i]/positionVector_cubed)};
        //std::vector<double> DipoleValues = { gConstant * (3.0 * mu1DotPosition * mu1DotPosition * positionVector[i] - mu1_dot_mu2 * positionVector[i]),
        //                                     gConstant * (3.0 * mu1DotPosition * mu1DotPosition * positionVector[i] - mu1_dot_mu2 * positionVector[i]),
        //                                     gConstant * (3.0 * mu1DotPosition * mu1DotPosition * positionVector[i] - mu1_dot_mu2 * positionVector[i])};

        //gConstant * (3.0 * mu1_dot_r12 * mu2_dot_r12 * positionVector[i] - mu1_dot_mu2 * positionVector[i]);
        //gConstant * ((3 * positionVector[2] * mu2DotPosition)/positionVector_fifth - mzTerms[i]/positionVector_cubed)
        totalDipoleTerms[0] += DipoleValues[0];
        totalDipoleTerms[1] += DipoleValues[1];
        totalDipoleTerms[2] += DipoleValues[2];
    }

    return totalDipoleTerms;
}

/*
std::vector<double> Numerical_Methods_Class::DipoleDipoleCoupling(double magneticMoment1, double magneticMoment2,
                                                                  int originSite, int targetSite) {
    double exchangeStiffness = 5.3e-17;
    double latticeConstant = sqrt(exchangeStiffness / _exchangeVec[originSite - 1]);

    std::vector<double> positionVector = {std::abs(targetSite - originSite)* latticeConstant, 0.0, 0.0};
    double positionVectorMagnitude = std::sqrt(positionVector[0] * positionVector[0] + positionVector[1] * positionVector[1] + positionVector[2] * positionVector[2]);

    double mu1_dot_r12 = magneticMoment1 * positionVector[0]; // + magneticMoment1[1] * positionVector[1] + magneticMoment1[2] * positionVector[2];
    double mu2_dot_r12 = magneticMoment2 * positionVector[0]; // + magneticMoment2[1] * positionVector[1] + magneticMoment2[2] * positionVector[2];
    double mu1_dot_mu2 = magneticMoment1 * magneticMoment2; // + magneticMoment1[1] * magneticMoment2[1] + magneticMoment1[2] * magneticMoment2[2];

    std::vector<double> hDipoleTerms(3);
    double gConstant = (_permFreeSpace / (4.0 * M_PI * std::pow(positionVectorMagnitude, 5)));
    for (int i = 0; i < 2; i++) {
        if (i == 1)
            continue;

        hDipoleTerms[i] =
                gConstant * (3.0 * mu1_dot_r12 * mu2_dot_r12 * positionVector[i] - mu1_dot_mu2 * positionVector[i]);
    }
    return hDipoleTerms;
}
*/

double Numerical_Methods_Class::EffectiveFieldX(const int& site, const double& mxLHS, const double& mxMID,
                                                const double& mxRHS, const double& dipoleTerm, const double& current_time) {
    // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
    double hx;

    if (_isFM) {
        if (site >= _drivingRegionLHS && site <= _drivingRegionRHS) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (_hasStaticDrive)
                hx = _exchangeVec[site - 1] * mxLHS + _exchangeVec[site] * mxRHS + dipoleTerm + _dynamicBiasField;
            else if (!_hasStaticDrive)
                hx = _exchangeVec[site - 1] * mxLHS + _exchangeVec[site] * mxRHS + dipoleTerm +
                      _dynamicBiasField * cos(_drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx = _exchangeVec[site - 1] * mxLHS + _exchangeVec[site] * mxRHS + dipoleTerm;
    } else if (!_isFM) {
        if (site >= _drivingRegionLHS && site <= _drivingRegionRHS) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if (_hasStaticDrive)
                hx = -1.0 * (_exchangeVec[site - 1] * mxLHS + _exchangeVec[site] * mxRHS + _dynamicBiasField);
            else if (!_hasStaticDrive)
                hx = -1.0 * (_exchangeVec[site - 1] * mxLHS + _exchangeVec[site] * mxRHS) + _dynamicBiasField * cos(_drivingAngFreq * current_time);
        } else
            // All spins along x which are not within the driving region
            hx = -1.0 * (_exchangeVec[site - 1] * mxLHS + _exchangeVec[site] * mxRHS);
    }

    return hx;
}
double Numerical_Methods_Class::EffectiveFieldY(const int& site, const double& myLHS, const double& myMID, const double& myRHS,
                                                const double &dipoleTerm) {
    // The effective field (H_eff) y-component acting upon a given magnetic moment (site), abbreviated to 'hy'
    double hy;

    if (_isFM) {
        hy = _exchangeVec[site - 1] * myLHS + _exchangeVec[site] * myRHS + dipoleTerm;
    } else if (!_isFM) {
        hy = -1.0 * (_exchangeVec[site - 1] * myLHS + _exchangeVec[site] * myRHS);
    }

    return hy;
}
double Numerical_Methods_Class::EffectiveFieldZ(const int& site, const double& mzLHS, const double& mzMID, const double& mzRHS,
                                                const double& dipoleTerm) {
    // The effective field (H_eff) z-component acting upon a given magnetic moment (site), abbreviated to 'hz'
    double hz;

    if (_isFM) {
        hz = _exchangeVec[site - 1] * mzLHS + _exchangeVec[site] * mzRHS + dipoleTerm + GV.GetStaticBiasField();
    } else if (!_isFM) {
        if (mzMID > 0)
            hz = GV.GetStaticBiasField() + _anisotropyField - (_exchangeVec[site - 1] * mzLHS + _exchangeVec[site] * mzRHS);
        else if (mzMID < 0)
            hz = GV.GetStaticBiasField() - _anisotropyField - (_exchangeVec[site - 1] * mzLHS + _exchangeVec[site] * mzRHS);
    }

    return hz;
}

double Numerical_Methods_Class::MagneticMomentX(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mxK; // The magnetic moment component along the x-direction for the first stage of RK2

    if (_useLLG) {
        // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK2.
        mxK = _gyroMagConst * (- (_gilbertVector[site] * hyMID * mxMID * myMID) + hyMID * mzMID - hzMID * (myMID + _gilbertVector[site] * mxMID * mzMID) + _gilbertVector[site] * hxMID * (pow(myMID,2) + pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mxK = -1.0 * _gyroMagConst * (myMID * hzMID - mzMID * hyMID);
    }

    return mxK;
}
double Numerical_Methods_Class::MagneticMomentY(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {
    double myK;

    if (_useLLG) {
        // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK2.
        myK = _gyroMagConst * (-(hxMID * mzMID) + hzMID * (mxMID - _gilbertVector[site] * myMID * mzMID) + _gilbertVector[site] * (hyMID * pow(mxMID,2) - hxMID * mxMID * myMID + hyMID * pow(mzMID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        myK = _gyroMagConst * (mxMID * hzMID - mzMID * hxMID);
    }

    return myK;
}
double Numerical_Methods_Class::MagneticMomentZ(const int& site, const double& mxMID, const double& myMID, const double& mzMID,
                                                const double& hxMID, const double& hyMID, const double& hzMID) {

    double mzK;

    if (_useLLG) {
        // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK2.
        mzK = _gyroMagConst * (hxMID * myMID + _gilbertVector[site] * hzMID * (pow(mxMID,2) + pow(myMID,2)) - _gilbertVector[site]*hxMID*mxMID*mzMID - hyMID * (mxMID + _gilbertVector[site] * myMID * mzMID));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mzK = -1.0 * _gyroMagConst * (mxMID * hyMID - myMID * hxMID);
    }

    return mzK;
}

void Numerical_Methods_Class::SolveRK2() {

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");

    if (_isFM) {
        InformUserOfCodeType("RK2 Midpoint (FM)");
        CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)");
    } else if (!_isFM) {
        InformUserOfCodeType("RK2 Midpoint (AFM)");
        CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
    }

    if (GV.GetEmailWhenCompleted()) {
        CreateMetadata();
    }

    progressbar bar(100);

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {

        if (_iterationEnd >= 100 && iteration % (_iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        TestShockwaveConditions(iteration);

        double t0 = _totalTime, t0HalfStep = _totalTime + _stepsizeHalf;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(GV.GetNumSpins() + 2, 0), my1(GV.GetNumSpins() + 2, 0), mz1(GV.GetNumSpins() + 2, 0);

        // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)
        // RK4 Stage 1. Takes initial conditions as inputs.

        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            double dipoleX, dipoleY, dipoleZ;
            if (_useDipolar) {
                std::vector<double> mxTermsForDipole = {_mx0[spinLHS], _mx0[site], _mx0[spinRHS]};
                std::vector<double> myTermsForDipole = {_my0[spinLHS], _my0[site], _my0[spinRHS]};
                std::vector<double> mzTermsForDipole = {_mz0[spinLHS], _mz0[site], _mz0[spinRHS]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipoleDipoleCoupling(mxTermsForDipole, myTermsForDipole,
                                                                       mzTermsForDipole, siteTermsForDipole);


                double dipoleX = dipoleTerms[0];
                double dipoleY = dipoleTerms[1];
                double dipoleZ = dipoleTerms[2];
            } else {
                dipoleX = 0;
                dipoleY = 0;
                dipoleZ = 0;
            }

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK0 = EffectiveFieldX(site, _mx0[spinLHS], _mx0[site], _mx0[spinRHS], dipoleX, t0);
            double hyK0 = EffectiveFieldY(site, _my0[spinLHS], _my0[site], _my0[spinRHS], dipoleY);
            double hzK0 = EffectiveFieldZ(site, _mz0[spinLHS], _mz0[site], _mz0[spinRHS], dipoleZ);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK1 = MagneticMomentX(site, _mx0[site], _my0[site], _mz0[site], hxK0, hyK0, hzK0);
            double myK1 = MagneticMomentY(site, _mx0[site], _my0[site], _mz0[site], hxK0, hyK0, hzK0);
            double mzK1 = MagneticMomentZ(site, _mx0[site], _my0[site], _mz0[site], hxK0, hyK0, hzK0);

            // Find (m0 + k1/2) for each site, which is used in the next stage.
            mx1[site] = _mx0[site] + _stepsizeHalf * mxK1;
            my1[site] = _my0[site] + _stepsizeHalf * myK1;
            mz1[site] = _mz0[site] + _stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        std::vector<double> mx2(GV.GetNumSpins() + 2, 0), my2(GV.GetNumSpins() + 2, 0), mz2(GV.GetNumSpins() + 2, 0);

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            double dipoleX, dipoleY, dipoleZ;
            if (_useDipolar) {
                std::vector<double> mxTermsForDipole = {mx1[spinLHS], mx1[site], mx1[spinRHS]};
                std::vector<double> myTermsForDipole = {my1[spinLHS], my1[site], my1[spinRHS]};
                std::vector<double> mzTermsForDipole = {mz1[spinLHS], mz1[site], mz1[spinRHS]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipoleDipoleCoupling(mxTermsForDipole, myTermsForDipole,
                                                                       mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            } else {
                dipoleX = 0;
                dipoleY = 0;
                dipoleZ = 0;
            }
            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK1 = EffectiveFieldX(site, mx1[spinLHS], mx1[site], mx1[spinRHS], dipoleX, t0);
            double hyK1 = EffectiveFieldY(site, my1[spinLHS], my1[site], my1[spinRHS], dipoleY);
            double hzK1 = EffectiveFieldZ(site, mz1[spinLHS], mz1[site], mz1[spinRHS], dipoleZ);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK2 = MagneticMomentX(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double myK2 = MagneticMomentY(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double mzK2 = MagneticMomentZ(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);

            mx2[site] = _mx0[site] + _stepsize * mxK2;
            my2[site] = _my0[site] + _stepsize * myK2;
            mz2[site] = _mz0[site] + _stepsize * mzK2;

            if (_shouldTrackMValues) {
                double mIterationNorm = sqrt(pow(mx2[site], 2) + pow(my2[site], 2) + pow(mz2[site], 2));
                if ((_largestMNorm) > (1.0 - mIterationNorm)) { _largestMNorm = (1.0 - mIterationNorm); }
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        _mx0.clear(); _my0.clear(); _mz0.clear();
        mx1.clear(); my1.clear(); mz1.clear();

        SaveDataToFile(mxRK2File, mx2, iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        _mx0 = mx2; _my0 = my2; _mz0 = mz2;

        if (iteration == _forceStopAtIteration)
            exit(0);

        _totalTime += _stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    if (GV.GetEmailWhenCompleted()) {
        CreateMetadata(true);
    }

    if (_shouldTrackMValues)
        std::cout << "\nMax norm. value of M is: " << _largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::SolveRK2Bilayer() {

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");

    if (_isFM) {
        InformUserOfCodeType("RK2 Midpoint (FM)");
        CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)");
    } else if (!_isFM) {
        InformUserOfCodeType("RK2 Midpoint (AFM)");
        CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
    }

    if (GV.GetEmailWhenCompleted()) {
        CreateMetadata();
    }

    progressbar bar(100);

    // Initialise mapping
    std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested1;

    // Assign name to nested-nested vector
    mValsNested1["nestedNestedVector1"] = initializeNestedNestedVector(1, true);

    // Assign name to nested-nested vector
    std::vector<std::vector<std::vector<double>>> m0Nest = mValsNested1["nestedNestedVector1"];

    // Invoke method to set initial magnetic moments. To call: mValsNest[layer][site][component]
    SetInitialMagneticMomentsMultilayer(m0Nest, 1, 0, 0 , 1);

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {

        if (_iterationEnd >= 100 && iteration % (_iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        TestShockwaveConditions(iteration);

        double t0 = _totalTime, t0HalfStep = _totalTime + _stepsizeHalf;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = mx0 + (h * k1 / 2) etc
        std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested2;
        mValsNested2["nestedNestedVector1"] = initializeNestedNestedVector(1, true);
        std::vector<std::vector<std::vector<double>>> m1Nest = mValsNested2["nestedNestedVector1"];
        SetInitialMagneticMomentsMultilayer(m1Nest, 1, 0, 0 , 0);

        // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)
        // RK2 Stage 1. Takes initial conditions as inputs.

        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            double dipoleX, dipoleY, dipoleZ;
            if (_useDipolar) {
                std::vector<double> mxTermsForDipole = {m0Nest[0][spinLHS][0], m0Nest[0][site][0], m0Nest[0][spinRHS][0]};
                std::vector<double> myTermsForDipole = {m0Nest[0][spinLHS][1], m0Nest[0][site][1], m0Nest[0][spinRHS][1]};
                std::vector<double> mzTermsForDipole = {m0Nest[0][spinLHS][2], m0Nest[0][site][2], m0Nest[0][spinRHS][2]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipoleDipoleCoupling(mxTermsForDipole, myTermsForDipole,
                                                                       mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            } else {
                dipoleX = 0;
                dipoleY = 0;
                dipoleZ = 0;
            }

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK0 = EffectiveFieldX(site, m0Nest[0][spinLHS][0], m0Nest[0][site][0], m0Nest[0][spinRHS][0], dipoleX, t0);
            double hyK0 = EffectiveFieldY(site, m0Nest[0][spinLHS][1], m0Nest[0][site][1], m0Nest[0][spinRHS][1], dipoleY);
            double hzK0 = EffectiveFieldZ(site, m0Nest[0][spinLHS][2], m0Nest[0][site][2], m0Nest[0][spinRHS][2], dipoleZ);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK1 = MagneticMomentX(site, m0Nest[0][site][0], m0Nest[0][site][1], m0Nest[0][site][2], hxK0, hyK0, hzK0);
            double myK1 = MagneticMomentY(site, m0Nest[0][site][0], m0Nest[0][site][1], m0Nest[0][site][2], hxK0, hyK0, hzK0);
            double mzK1 = MagneticMomentZ(site, m0Nest[0][site][0], m0Nest[0][site][1], m0Nest[0][site][2], hxK0, hyK0, hzK0);

            // Find (m0 + k1/2) for each site, which is used in the next stage.
            m1Nest[0][site][0] =  m0Nest[0][site][0] + _stepsizeHalf * mxK1;
            m1Nest[0][site][1] =  m0Nest[0][site][1] + _stepsizeHalf * myK1;
            m1Nest[0][site][2] =  m0Nest[0][site][2] + _stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested3;
        mValsNested3["nestedNestedVector1"] = initializeNestedNestedVector(1, true);
        std::vector<std::vector<std::vector<double>>> m2Nest = mValsNested3["nestedNestedVector1"];
        SetInitialMagneticMomentsMultilayer(m2Nest, 1, 0, 0 , 0);

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            double dipoleX, dipoleY, dipoleZ;
            if (_useDipolar) {
                std::vector<double> mxTermsForDipole = {m1Nest[0][spinLHS][0], m1Nest[0][site][0], m1Nest[0][spinRHS][0]};
                std::vector<double> myTermsForDipole = {m1Nest[0][spinLHS][1], m1Nest[0][site][1], m1Nest[0][spinRHS][1]};
                std::vector<double> mzTermsForDipole = {m1Nest[0][spinLHS][2], m1Nest[0][site][2], m1Nest[0][spinRHS][2]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipoleDipoleCoupling(mxTermsForDipole, myTermsForDipole,
                                                                       mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            } else {
                dipoleX = 0;
                dipoleY = 0;
                dipoleZ = 0;
            }
            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK1 = EffectiveFieldX(site, m1Nest[0][spinLHS][0], m1Nest[0][site][0], m1Nest[0][spinRHS][0], dipoleX, t0);
            double hyK1 = EffectiveFieldY(site, m1Nest[0][spinLHS][1], m1Nest[0][site][1], m1Nest[0][spinRHS][1], dipoleY);
            double hzK1 = EffectiveFieldZ(site, m1Nest[0][spinLHS][2], m1Nest[0][site][2], m1Nest[0][spinRHS][2], dipoleZ);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK2 = MagneticMomentX(site, m1Nest[0][site][0], m1Nest[0][site][1], m1Nest[0][site][2], hxK1, hyK1, hzK1);
            double myK2 = MagneticMomentY(site, m1Nest[0][site][0], m1Nest[0][site][1], m1Nest[0][site][2], hxK1, hyK1, hzK1);
            double mzK2 = MagneticMomentZ(site, m1Nest[0][site][0], m1Nest[0][site][1], m1Nest[0][site][2], hxK1, hyK1, hzK1);

            m2Nest[0][site][0] = m0Nest[0][site][0] + _stepsize * mxK2;
            m2Nest[0][site][1] = m0Nest[0][site][1] + _stepsize * myK2;
            m2Nest[0][site][2] = m0Nest[0][site][2] + _stepsize * mzK2;

            if (_shouldTrackMValues) {
                double mIterationNorm = sqrt(pow(m2Nest[0][site][0], 2) + pow(m2Nest[0][site][1], 2) + pow(m2Nest[0][site][2], 2));
                if ((_largestMNorm) > (1.0 - mIterationNorm)) { _largestMNorm = (1.0 - mIterationNorm); }
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        m0Nest.shrink_to_fit();
        m1Nest.shrink_to_fit();

        SaveDataToFileMultilayer(mxRK2File, m2Nest[0], iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        m0Nest = m2Nest;

        if (iteration == _forceStopAtIteration)
            exit(0);

        _totalTime += _stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    if (GV.GetEmailWhenCompleted()) {
        CreateMetadata(true);
    }

    if (_shouldTrackMValues)
        std::cout << "\nMax norm. value of M is: " << _largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::SolveRK2BilayerTest() {

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");

    if (_isFM) {
        InformUserOfCodeType("RK2 Midpoint (FM)");
        CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)");
    } else if (!_isFM) {
        InformUserOfCodeType("RK2 Midpoint (AFM)");
        CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
    }

    if (GV.GetEmailWhenCompleted()) {
        CreateMetadata();
    }

    progressbar bar(100);

    // Initialise mapping
    std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested1;

    // Assign name to nested-nested vector
    mValsNested1["nestedNestedVector1"] = initializeNestedNestedVector(1, true);

    // Assign name to nested-nested vector
    std::vector<std::vector<std::vector<double>>> m0Nest = mValsNested1["nestedNestedVector1"];

    // Invoke method to set initial magnetic moments. To call: mValsNest[layer][site][component]
    SetInitialMagneticMomentsMultilayer(m0Nest, 1, 0, 0 , 1);

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {

        if (_iterationEnd >= 100 && iteration % (_iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        TestShockwaveConditions(iteration);

        double t0 = _totalTime, t0HalfStep = _totalTime + _stepsizeHalf;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = mx0 + (h * k1 / 2) etc
        std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested2;
        mValsNested2["nestedNestedVector1"] = initializeNestedNestedVector(1, true);
        std::vector<std::vector<std::vector<double>>> m1Nest = mValsNested2["nestedNestedVector1"];
        SetInitialMagneticMomentsMultilayer(m1Nest, 1, 0, 0 , 0);

        // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)
        // RK2 Stage 1. Takes initial conditions as inputs.

        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            double mxLHS = m0Nest[0][spinLHS][0], mxMID = m0Nest[0][site][0], mxRHS = m0Nest[0][spinRHS][0];
            double myLHS = m0Nest[0][spinLHS][1], myMID = m0Nest[0][site][1], myRHS = m0Nest[0][spinRHS][1];
            double mzLHS = m0Nest[0][spinLHS][2], mzMID = m0Nest[0][site][2], mzRHS = m0Nest[0][spinRHS][2];

            double dipoleX, dipoleY, dipoleZ;
            if (_useDipolar) {
                std::vector<double> mxTermsForDipole = {m0Nest[0][spinLHS][0], m0Nest[0][site][0], m0Nest[0][spinRHS][0]};
                std::vector<double> myTermsForDipole = {m0Nest[0][spinLHS][1], m0Nest[0][site][1], m0Nest[0][spinRHS][1]};
                std::vector<double> mzTermsForDipole = {m0Nest[0][spinLHS][2], m0Nest[0][site][2], m0Nest[0][spinRHS][2]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipoleDipoleCoupling(mxTermsForDipole, myTermsForDipole,
                                                                       mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            } else {
                dipoleX = 0;
                dipoleY = 0;
                dipoleZ = 0;
            }

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK0 = EffectiveFieldX(site, mxLHS, mxMID, mxRHS, dipoleX, t0);
            double hyK0 = EffectiveFieldY(site, myLHS, myMID, myRHS, dipoleY);
            double hzK0 = EffectiveFieldZ(site, mzLHS, mzMID, mzRHS, dipoleZ);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK1 = MagneticMomentX(site, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
            double myK1 = MagneticMomentY(site, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
            double mzK1 = MagneticMomentZ(site, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);

            // Find (m0 + k1/2) for each site, which is used in the next stage.
            m1Nest[0][site][0] =  m0Nest[0][site][0] + _stepsizeHalf * mxK1;
            m1Nest[0][site][1] =  m0Nest[0][site][1] + _stepsizeHalf * myK1;
            m1Nest[0][site][2] =  m0Nest[0][site][2] + _stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        std::map<std::string, std::vector<std::vector<std::vector<double>>>> mValsNested3;
        mValsNested3["nestedNestedVector1"] = initializeNestedNestedVector(1, true);
        std::vector<std::vector<std::vector<double>>> m2Nest = mValsNested3["nestedNestedVector1"];
        SetInitialMagneticMomentsMultilayer(m2Nest, 1, 0, 0 , 0);

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            double mxLHS = m1Nest[0][spinLHS][0], mxMID = m1Nest[0][site][0], mxRHS = m1Nest[0][spinRHS][0];
            double myLHS = m1Nest[0][spinLHS][1], myMID = m1Nest[0][site][1], myRHS = m1Nest[0][spinRHS][1];
            double mzLHS = m1Nest[0][spinLHS][2], mzMID = m1Nest[0][site][2], mzRHS = m1Nest[0][spinRHS][2];

            double dipoleX, dipoleY, dipoleZ;
            if (_useDipolar) {
                std::vector<double> mxTermsForDipole = {m1Nest[0][spinLHS][0], m1Nest[0][site][0], m1Nest[0][spinRHS][0]};
                std::vector<double> myTermsForDipole = {m1Nest[0][spinLHS][1], m1Nest[0][site][1], m1Nest[0][spinRHS][1]};
                std::vector<double> mzTermsForDipole = {m1Nest[0][spinLHS][2], m1Nest[0][site][2], m1Nest[0][spinRHS][2]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipoleDipoleCoupling(mxTermsForDipole, myTermsForDipole,
                                                                       mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            } else {
                dipoleX = 0;
                dipoleY = 0;
                dipoleZ = 0;
            }
            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK1 = EffectiveFieldX(site, mxLHS, mxMID, mxRHS, dipoleX, t0);
            double hyK1 = EffectiveFieldY(site, myLHS, myMID, myRHS, dipoleY);
            double hzK1 = EffectiveFieldZ(site, mzLHS, mzMID, mzRHS, dipoleZ);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK2 = MagneticMomentX(site, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
            double myK2 = MagneticMomentY(site, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
            double mzK2 = MagneticMomentZ(site, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);

            m2Nest[0][site][0] = m0Nest[0][site][0] + _stepsize * mxK2;
            m2Nest[0][site][1] = m0Nest[0][site][1] + _stepsize * myK2;
            m2Nest[0][site][2] = m0Nest[0][site][2] + _stepsize * mzK2;

            if (_shouldTrackMValues) {
                double mIterationNorm = sqrt(pow(m2Nest[0][site][0], 2) + pow(m2Nest[0][site][1], 2) + pow(m2Nest[0][site][2], 2));
                if ((_largestMNorm) > (1.0 - mIterationNorm)) { _largestMNorm = (1.0 - mIterationNorm); }
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */

        for(auto& v : m0Nest) {
            v.clear();
            v.shrink_to_fit();
        }
        for(auto& v : m1Nest) {
            v.clear();
            v.shrink_to_fit();
        }
        m0Nest.clear(); m0Nest.shrink_to_fit();
        m1Nest.clear(); m1Nest.shrink_to_fit();

        SaveDataToFileMultilayer(mxRK2File, m2Nest[0], iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        m0Nest = m2Nest;

        if (iteration == _forceStopAtIteration)
            exit(0);

        _totalTime += _stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    if (GV.GetEmailWhenCompleted()) {
        CreateMetadata(true);
    }

    if (_shouldTrackMValues)
        std::cout << "\nMax norm. value of M is: " << _largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata) {
    /**
     * Write all non-data information to the output file.
     */
    if (is_metadata) {
        outputFileName << "Key Data\n\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using LLG: [" << _useLLG << "]\t\t\t\tUsing Shockwave: [" << _hasShockwave << "]\t\tDrive from LHS: [" << _lhsDrive <<
                       "]\nNumerical Method Used: [" << methodUsed << "]\t\tHas Static Drive: [" << _hasStaticDrive << "]\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0): " << GV.GetStaticBiasField() << " T\t\t\t" << "Dynamic Bias Field (H_D1): " << _dynamicBiasField << " T\n" <<
                          "Dynamic Bias Field Scale Factor: " << _shockwaveInitialStrength << "\t\t" << "Second Dynamic Bias Field (H_D2): " << _shockwaveMax << " T\n" <<
                          "Driving Frequency (f): " << _drivingFreq << "Hz\t\t""Driving Region Start Site: " << _drivingRegionLHS - _numSpinsDamped << "\n" <<
                          "Driving Region End Site: " << _drivingRegionRHS - _numSpinsDamped << " \t\t\t" << "Driving Region Width: " << _drivingRegionWidth << " \n" <<
                          "Max. Sim. Time: " << _maxSimTime << " s\t\t\t\t" << "Min. Exchange Val (J): " << GV.GetExchangeMinVal()  << " T\n" <<
                          "Max. Exchange Val (J): " << GV.GetExchangeMaxVal() << " T\t\t\t" << "Max. Iterations: " << _iterationEnd << "\n" <<
                          "No. DataPoints: " << _numberOfDataPoints << " \t\t\t\t" << "No. Spins in Chain: " << _numSpinsInChain << "\n" <<
                          "No. Damped Spins: " << _numSpinsDamped << "per side\t\t\t" << "No. Total Spins: " << GV.GetNumSpins() << " \n" <<
                          "Stepsize (h): " << _stepsize << "\t\t\t\t" << "Gilbert Damping Factor: " << _gilbertConst << "\n" <<
                          "Gyromagnetic Ratio (2Pi*Y): " << _gyroMagConst << "\t\t""Shockwave Gradient Time: " << _iterStartShock << "s\n" <<
                          "Shockwave Application Time: " << _shockwaveGradientTime * _stepsize << "s\n" <<
                          std::endl;

        return;
    }
    else {

        outputFileName << "Key Data\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using LLG," << _useLLG << ",Using Shockwave," << _hasShockwave << ",Drive from LHS," << _lhsDrive <<
                       ",Numerical Method Used," << methodUsed << ",Has Static Drive," << _hasStaticDrive << "\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                          "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site, Driving Region Width,"
                          "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                          "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, Stepsize (h),Gilbert Damping Factor, Gyromagnetic Ratio (2Pi*Y),"
                          "Shockwave Gradient Time [s], Shockwave Application Time [s]"
                          "\n";

        outputFileName << GV.GetStaticBiasField() << ", " << _dynamicBiasField << ", " << _shockwaveInitialStrength << ", " << _shockwaveMax << ", "
                       << _drivingFreq << ", " << _drivingRegionLHS - _numSpinsDamped << ", " << _drivingRegionRHS - _numSpinsDamped << ", " << _drivingRegionWidth << ", "
                       << _maxSimTime << ", " << GV.GetExchangeMinVal() << ", " << GV.GetExchangeMaxVal() << ", " << _iterationEnd << ", " << _numberOfDataPoints << ", "
                       << _numSpinsInChain << ", " << _numSpinsDamped << ", " << GV.GetNumSpins() << ", " << _stepsize << ", " << _gilbertConst << ", " << _gyroMagConst << ", "
                       << _iterStartShock << ", " << _shockwaveGradientTime * _stepsize
                       << "\n";

        outputFileName << "\n";
    }

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments );
    outputFileName << "Note(s):," << notesComments << "\n"; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n";

    outputFileName << "\n";

    CreateColumnHeaders(outputFileName);

    std::cout << "\n";
}
void Numerical_Methods_Class::CreateColumnHeaders(std::ofstream &outputFileName) {
    /**
     * Creates the column headers for each spin site simulated. This code can change often, so compartmentalising it in
     * a separate function is necessary to reduce bugs.
     */
    if (_printAllData or _printFixedLines) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for (int i = 1; i <= GV.GetNumSpins(); i++) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if (_printFixedSites) {

        outputFileName << "Time";
        for (int & fixed_out_val : _fixed_output_sites)
            outputFileName << "," << fixed_out_val;
        outputFileName << std::endl;

        //outputFileName << "Time" << ", "
        //               << static_cast<int>(14000) << ","
        //               << static_cast<int>(16000) << ","
        //               << static_cast<int>(18000) << ","
        //               << static_cast<int>(20000) << std::endl;

    }
}
void Numerical_Methods_Class::InformUserOfCodeType(const std::string& nameNumericalMethod) {
    /**
     * Informs the user of the code type they are running, including: solver type; special modules.
     */
    if (_useLLG)
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (LLG) code";
    else
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (Torque) code";

    if (_hasShockwave)
        std::cout << " with shockwave module.\n";
    else
        std::cout << ".\n";

}
void Numerical_Methods_Class::PrintVector(std::vector<double> &vectorToPrint, bool shouldExitAfterPrint) {

    std::cout << "\n\n";

    int count = 0;
    for (double i: vectorToPrint) {
        if (++count % 10 == 0)
            std::cout << std::setw(12) << i << std::endl;
        else
            std::cout << std::setw(12) << i << ", ";
    }

    if (shouldExitAfterPrint)
        exit(0);
}
void Numerical_Methods_Class::PrintNestedNestedVector(std::vector<std::vector<std::vector<double>>> nestedNestedVector){
        // Print the contents of nestedNestedVector1
    std::cout << "{";
    for (size_t i = 0; i < nestedNestedVector.size(); ++i) {
        std::cout << "{";
        for (size_t j = 0; j < nestedNestedVector[i].size(); ++j) {
            std::cout << "{";
            for (size_t k = 0; k < nestedNestedVector[i][j].size(); ++k) {
                std::cout << nestedNestedVector[i][j][k];
                if (k < nestedNestedVector[i][j].size() - 1) {
                    std::cout << ",";
                }
            }
            std::cout << "}";
            if (j < nestedNestedVector[i].size() - 1) {
                std::cout << ",";
            }
        }
        std::cout << "}";
        if (i < nestedNestedVector.size() - 1) {
            std::cout << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}
void Numerical_Methods_Class::SaveDataToFile(std::ofstream &outputFileName, std::vector<double> &arrayToWrite, int &iteration) {
    std::cout.precision(6);
    std::cout << std::scientific;

    if (iteration % (_iterationEnd / _numberOfDataPoints) == 0) {
        if (_printFixedLines) {
            for (int i = 0; i <= GV.GetNumSpins(); i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * _stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if (_printFixedSites) {
            /*outputFileName << (iteration * _stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * _stepsize);
            for (int & fixed_out_val : _fixed_output_sites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (_printAllData) {
        for (int i = 0; i <= GV.GetNumSpins(); i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * _stepsize) << ","; // Print current time
            else if (i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (_printFixedLines) {
        // iteration >= static_cast<int>(_iterationEnd / 2.0) &&
        if (iteration % (_iterationEnd / _numberOfDataPoints) == 0) {
            //if (iteration == _iterationEnd) {
            for (int i = 0; i <= GV.GetNumSpins(); i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * _stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;
        }
    } else {
        if (_printAllData) {
            for (int i = 0; i <= GV.GetNumSpins(); i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * _stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (_iterationEnd / _numberOfDataPoints) == 0) {
                if (_printFixedSites) {

                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[_drivingRegionLHS] << ","
                                   << arrayToWrite[static_cast<int>(_drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[_drivingRegionRHS] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[GV.GetNumSpins()] << std::endl;

                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[_drivingRegionLHS] << ","
                                   << arrayToWrite[static_cast<int>(_drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[_drivingRegionRHS] << ","
                                   << arrayToWrite[static_cast<int>(GV.GetNumSpins() / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(GV.GetNumSpins() / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(GV.GetNumSpins() / 4.0)] << ","
                                   << arrayToWrite[GV.GetNumSpins()] << std::endl;
                }
            }
        }
    } */
}
void Numerical_Methods_Class::SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::vector<double> arrayToWrite;
    // Extract the first element from each nested vector
    for (const auto& innerVector : nestedArrayToWrite) {
        if (!innerVector.empty()) {
            arrayToWrite.push_back(innerVector[0]);
        }
    }

    if (iteration % (_iterationEnd / _numberOfDataPoints) == 0) {
        if (_printFixedLines) {
            for (int i = 0; i <= GV.GetNumSpins(); i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * _stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if (_printFixedSites) {
            /*outputFileName << (iteration * _stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * _stepsize);
            for (int & fixed_out_val : _fixed_output_sites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (_printAllData) {
        for (int i = 0; i <= GV.GetNumSpins(); i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * _stepsize) << ","; // Print current time
            else if (i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (_printFixedLines) {
        // iteration >= static_cast<int>(_iterationEnd / 2.0) &&
        if (iteration % (_iterationEnd / _numberOfDataPoints) == 0) {
            //if (iteration == _iterationEnd) {
            for (int i = 0; i <= GV.GetNumSpins(); i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * _stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;
        }
    } else {
        if (_printAllData) {
            for (int i = 0; i <= GV.GetNumSpins(); i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * _stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (_iterationEnd / _numberOfDataPoints) == 0) {
                if (_printFixedSites) {

                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[_drivingRegionLHS] << ","
                                   << arrayToWrite[static_cast<int>(_drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[_drivingRegionRHS] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[GV.GetNumSpins()] << std::endl;

                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[_drivingRegionLHS] << ","
                                   << arrayToWrite[static_cast<int>(_drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[_drivingRegionRHS] << ","
                                   << arrayToWrite[static_cast<int>(GV.GetNumSpins() / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(GV.GetNumSpins() / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(GV.GetNumSpins() / 4.0)] << ","
                                   << arrayToWrite[GV.GetNumSpins()] << std::endl;
                }
            }
        }
    } */
}
void Numerical_Methods_Class::TestShockwaveConditions(double iteration) {

    if (_shouldDriveCease) {
        // and (_isShockwaveOn and _isShockwaveAtMax)) {
        if (_isShockwaveOn and not _isShockwaveAtMax) {
            std::cout << "Shock not at maximum when cut-off" << std::endl;
        }

        if (iteration >= _iterationEnd * _iterEndShock) {
            // Shockwave begins once simulation is a certain % complete
            _hasShockwave = false;
            _isShockwaveOn = false;
            _dynamicBiasField = 0;
        }

        return;

    }

    // If method is triggered, then the applied biasFieldDriving is increased by the scale factor _shockwaveScaling
    if (_hasShockwave and not _isShockwaveOn)
    {
        if (iteration >= _iterationEnd * _iterStartShock)
        {
            // Shockwave begins once simulation is a certain % complete
            _isShockwaveOn = true;
            _dynamicBiasField = _shockwaveInitialStrength;
        }

        return;
    }

    if (_isShockwaveOn and not _isShockwaveAtMax)
    {
        _dynamicBiasField += _shockwaveStepsize;

        if (_dynamicBiasField >= _shockwaveMax)
        {
            _dynamicBiasField = _shockwaveMax;
            _isShockwaveAtMax = true;

        }
        return;

    }

}
void Numerical_Methods_Class::CreateMetadata(bool print_end_time) {

    std::string file_name = "simulation_metadata.txt";

    if (print_end_time) {
        std::ofstream metadata_end;
        metadata_end.open(GV.GetFilePath() + file_name, std::ios_base::app); // append instead of overwrite
        auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        metadata_end << "Finished at:\t" << std::put_time(localtime(&end), "%F %H-%M-%S") << std::endl;
        metadata_end.close();
    }
    else {
        std::ofstream metadata_start(GV.GetFilePath() + file_name);
        CreateFileHeader(metadata_start, "NM 2", true);
        auto start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        metadata_start << "Started at:\t" << std::put_time(localtime(&start), "%F %H-%M-%S") << std::endl;
        metadata_start.close();
    }
}
