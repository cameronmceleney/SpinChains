#include "Numerical_Methods_Class.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>

void Numerical_Methods_Class::NMSetup() {

    // ###################### Flags ######################
    _hasShockwave = false;
    _hasStaticDrive = false;
    _isFM = GV.GetIsFerromagnetic();
    _lhsDrive = true;
    _centralDrive = false;
    _dualDrive = false;
    _shouldTrackMValues = true;
    _useLLG = true;

    // ###################### Core Parameters ######################
    _drivingFreq = 12.5  * 1e9;
    _dynamicBiasField = 3e-4;
    _forceStopAtIteration = -1;
    _gyroMagConst = GV.GetGyromagneticConstant();
    _iterationEnd = static_cast<int>(7e5);  // 1e8
    _stepsize = 1e-15;  // 1e-17

    // ###################### Shockwave Parameters ######################
    _iterStartShock = 0.0;
    _shockwaveScaling = 1;
    _shockwaveGradientTime = 5e3;
    _shockwaveInitialStrength = 0;  // Set equal to _dynamicBiasField if NOT starting at time=0
    _shockwaveMax = 3e-3;

    // ###################### Data Output Parameters ######################
    _numberOfDataPoints = 1e2;  // 1e7
    _fixed_output_sites = {14000, 16000, 18000, 20000};
    _printFixedSites = false;
    _printFixedLines = true;
    _saveAllSpins = false;

    // ###################### Damping Factors ######################
    _gilbertConst  = 1e-6;
    _gilbertLower = 1e-6;
    _gilbertUpper = 1e0;

    // ###################### SpinChain Length Parameters ######################
    _drivingRegionWidth = 10;
    _numSpinsDamped = 0;

    // ###################### Computations based upon other inputs ######################
    _drivingAngFreq = 2 * M_PI * _drivingFreq;
    _maxSimTime = _stepsize * _iterationEnd;
    _numSpinsInChain = GV.GetNumSpins();
    _numberOfSpinPairs = _numSpinsInChain - 1;
    GV.SetNumSpins(_numSpinsInChain + 2 * _numSpinsDamped);
    _stepsizeHalf = _stepsize / 2.0;

    if (_isFM)
        _anisotropyField = 0;
    else if (!_isFM)
        _anisotropyField = GV.GetAnisotropyField();

    // ###################### Core Method Invocations ######################
    // Order is intentional and must be maintained!
    SetShockwaveConditions();
    SetDampingRegion();
    SetDrivingRegion();
    SetExchangeVector();
    SetInitialMagneticMoments();
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
    if ((_lhsDrive && _centralDrive) || (_lhsDrive && _dualDrive) || (_centralDrive && _dualDrive)) {
        std::cout << "Two (or more) conflicting driving region booleans were TRUE"
                  << "\n_lhsDrive: " << _lhsDrive << "\n_centralDrive: " << _centralDrive << "\n_dualDrive: " << _dualDrive
                  << "\n\nExiting...";
        exit(1);
    }

    if (_centralDrive) {
        _drivingRegionLHS = (_numSpinsInChain/2) +_numSpinsDamped - (_drivingRegionWidth / 2);
        _drivingRegionRHS = (_numSpinsInChain/2) +_numSpinsDamped + (_drivingRegionWidth / 2);
        return;
    }

    if (_dualDrive) {
        // TODO: The IF statement for _dualdrive still needs to be written
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

void Numerical_Methods_Class::RK2OriginalFM() {

    SpinChainEigenSolverClass printtest;

    // Notifies the user of what code they are running
    std::cout << "\nYou are running the RK2 Spinchains code." << std::endl;

    // Create files to save the data. All files will have (namefile) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath()+"rk2_mx_"+GV.GetFileNameBase()+".csv");
    std::ofstream myRK2File(GV.GetFilePath()+"rk2_my_"+GV.GetFileNameBase()+".csv");
    std::ofstream mzRK2File(GV.GetFilePath()+"rk2_mz_"+GV.GetFileNameBase()+".csv");

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or halfstep (h/2)*/
    for (long iteration = static_cast<long>(_iterationStart); iteration <= (long) _iterationEnd; iteration++) {

        _totalTime += _stepsize;
        double t0 = _totalTime; // The initial time of the iteration, and the time at the first stage; the start of the interval and the first step of RK2. Often called 't0' in literature
        double t0HalfStep = _totalTime + _stepsizeHalf; // The time at the midpoint of the interval; where the second stage of RK2 is. Often called 't0*h/2' in literature

        // The estimate of the slope for the x-axis magnetic moment component at the midpoint
        std::vector<double> mxEstMid(GV.GetNumSpins()+2, 0);
        // The estimate of the slope for the y-axis magnetic moment component at the midpoint
        std::vector<double> myEstMid(GV.GetNumSpins()+2, 0);
        // The estimate of the slope for the z-axis magnetic moment component at the midpoint
        std::vector<double> mzEstMid(GV.GetNumSpins()+2, 0);

        // Loop the 0th and final spins as they will always be zero-valued
        for (int spin = 1; spin <= GV.GetNumSpins()+1; spin++) {
            /* The first stage is based upon finding the value of the slope at the beginning of the interval (k1). This
             * stage takes the start conditions as an input, and substitutes them into the LLG equation. */
            int LHS_spin = spin - 1, RHS_spin = spin + 1;

            // The mx components for the first stage for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mx1 = _mx0[spin], mx1LHS = _mx0[LHS_spin], mx1RHS = _mx0[RHS_spin];
            // The my components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double my1 = _my0[spin], my1LHS = _my0[LHS_spin], my1RHS = _my0[RHS_spin];
            // The mz components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mz1 = _mz0[spin], mz1LHS = _mz0[LHS_spin], mz1RHS = _mz0[RHS_spin];

            double mx1K1, my1K1, mz1K1; // These are the estimations of the slopes at the beginning of the interval for each magnetic moment component
            double HeffX1K1, HeffY1K1, HeffZ1K1; // The effective field component acting upon each spin

            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                HeffX1K1 = _exchangeVec[LHS_spin] * mx1LHS + _exchangeVec[spin] * mx1RHS + _dynamicBiasField*cos(_drivingAngFreq * t0);
            } else {
                // The else statement includes all spins along x which are not within the driving region
                HeffX1K1 = _exchangeVec[LHS_spin] * mx1LHS + _exchangeVec[spin] * mx1RHS;
            }

            // No changes are made to the effective field in the y-direction
            HeffY1K1 = _exchangeVec[LHS_spin] * my1LHS + _exchangeVec[spin] * my1RHS;
            // The bias field is applied in the z-direction and so it contributes to the effective field in the z-direction
            HeffZ1K1 = _exchangeVec[LHS_spin] * mz1LHS + _exchangeVec[spin] * mz1RHS + GV.GetStaticBiasField();

            /* The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the
             * first stage of RK2.*/
            mx1K1 = -1 * _gyroMagConst * (my1 * HeffZ1K1 - mz1 * HeffY1K1);
            my1K1 = +1 * _gyroMagConst * (mx1 * HeffZ1K1 - mz1 * HeffX1K1);
            mz1K1 = -1 * _gyroMagConst * (mx1 * HeffY1K1 - my1 * HeffX1K1);

            mxEstMid[spin] = mx1 + mx1K1*_stepsizeHalf;
            myEstMid[spin] = my1 + my1K1*_stepsizeHalf;
            mzEstMid[spin] = mz1 + mz1K1*_stepsizeHalf;
        }

        // The estimate of the mx value for the next iteration of iteration calculated using the RK2 Midpoint rule
        std::vector<double> mxNextVal(GV.GetNumSpins()+2,0);
        // The estimate of the my value for the next iteration of iteration
        std::vector<double> myNextVal(GV.GetNumSpins()+2,0);
        // The estimate of the mz value for the next iteration of iteration
        std::vector<double> mzNextVal(GV.GetNumSpins()+2,0);

        for (int spin = 1; spin <= GV.GetNumSpins()+1; spin++) {
            /* The second stage uses the previously found k1 value, as well as the initial conditions, to determine the
             * value of the slope (k2) at the midpoint. In RK2, the values of k1 and k2 can then be jointly used to
             * estimate the next point of the function through a weighted average of k1 & k2.
             *
             * In this loop the definitions of variables follow a similar format to Stage1.*/

            int LHS_spin = spin - 1, RHS_spin = spin + 1;
            double mx2 = mxEstMid[spin], mx2LHS = mxEstMid[LHS_spin], mx2RHS = mxEstMid[RHS_spin];
            double my2 = myEstMid[spin], my2LHS = myEstMid[LHS_spin], my2RHS = myEstMid[RHS_spin];
            double mz2 = mzEstMid[spin], mz2LHS = mzEstMid[LHS_spin], mz2RHS = mzEstMid[RHS_spin];

            double mx2K2, my2K2, mz2K2;
            double HeffX2K2, HeffY2K2, HeffZ2K2;

            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // Driving region must be consistently applied at every stage of the RK2 method
                HeffX2K2 = _exchangeVec[LHS_spin] * mx2LHS + _exchangeVec[spin] * mx2RHS + _dynamicBiasField*cos(_drivingAngFreq * t0HalfStep);
            } else {
                HeffX2K2 = _exchangeVec[LHS_spin] * mx2LHS + _exchangeVec[spin] * mx2RHS;
            }

            HeffY2K2 = _exchangeVec[LHS_spin] * my2LHS + _exchangeVec[spin] * my2RHS;
            HeffZ2K2 = _exchangeVec[LHS_spin] * mz2LHS + _exchangeVec[spin] * mz2RHS + GV.GetStaticBiasField();

            mx2K2 = -1 * _gyroMagConst * ( my2*HeffZ2K2 - mz2*HeffY2K2 );
            my2K2 = +1 * _gyroMagConst * ( mx2*HeffZ2K2 - mz2*HeffX2K2 );
            mz2K2 = -1 * _gyroMagConst * ( mx2*HeffY2K2 - my2*HeffX2K2 );

            mxNextVal[spin] = _mx0[spin] + mx2K2*_stepsize;
            myNextVal[spin] = _my0[spin] + my2K2*_stepsize;
            mzNextVal[spin] = _mz0[spin] + mz2K2*_stepsize;

        } // Final line of the RK2 solvers for this iteration. Everything below here is part of the class function, but not the internal RK2 stage loops

        // Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear these between
        // loop iterations sometimes led to incorrect values cropping up
        _mx0.clear();
        _my0.clear();
        _mz0.clear();
        mxEstMid.clear();
        myEstMid.clear();
        mzEstMid.clear();

        /* Output function to write magnetic moment components to the terminal and/or files. Modulus component of IF
         * statement (default: 0.01 indicates how often the writing should occur. A value of 0.01 would mean writing
         * should occur every 1% of progress through the simulation*/
        if ( iteration % int(_iterationEnd*0.01) == 0 ) {
            std::cout << "Reporting Point: " << iteration << " iterations." << std::endl;
            /*  Code this is useful for debugging
             *  std::cout << "Numspins: " << GV.GetNumSpins() << "; const Jmin: " << GV.GetExchangeMinVal() << "; const Jmax: " << GV.GetExchangeMaxVal() << std::endl;
             *  std::cout << "RegionLHS: " << _drivingRegionLHS << "; RegionWidth: " << _drivingRegionWidth << "; RegionRHS: " << _drivingRegionRHS << std::endl;
             *  std::cout << "StepSize: " << _stepsize << "; HalfStepSize: " << _stepsizeHalf << "; TotalTime: " << _totalTime << "\n\n" << std::endl;
             *  printtest.PrintVector(mxNextVal); */

            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            for (int j = 1; j < GV.GetNumSpins() + 1; j++) {
                mxRK2File << mxNextVal[j] << ",";
                myRK2File << myNextVal[j] << ",";
                mzRK2File << mzNextVal[j] << ",";

                // Ensures that the final line doesn't contain a comma
                if (j == GV.GetNumSpins()) {
                    mxRK2File << mxNextVal[j] << std::flush;
                    myRK2File << myNextVal[j] << std::flush;
                    mzRK2File << mzNextVal[j] << std::flush;
                }
            }

            // Housekeeping
            mxRK2File << std::endl;
            myRK2File << std::endl;
            mzRK2File << std::endl;
        }

        /* Sets the final value of the current iteration of the loop (y_(n+1) in textbook's notation) to be the starting
         * value of the next iteration (y_n) */
        _mx0 = mxNextVal;
        _my0 = myNextVal;
        _mz0 = mzNextVal;
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();
    myRK2File.close();
    mzRK2File.close();

    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "Finished RK2 with: stepSize = " << _stepsize << "; itermax = " << _iterationEnd << "; filename = " << GV.GetFileNameBase() <<  std::endl;
}

void Numerical_Methods_Class::RK2MidpointFM() {

    progressbar bar(100);

    InformUserOfCodeType("RK2 Midpoint (FM)");

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)");

    if (GV.GetEmailWhenCompleted()) {
        CreateMetadata();
    }

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {

        if (_iterationEnd >= 100 && iteration % (_iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        TestShockwaveConditions(iteration);

        double t0 = _totalTime, t0HalfStep = _totalTime + _stepsizeHalf;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(GV.GetNumSpins() + 2, 0), my1(GV.GetNumSpins() + 2, 0), mz1(GV.GetNumSpins() + 2, 0);

        // Excludes the 0th and last spins as they will always be zero-valued (end, pinned spins)
        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            // RK2 Stage 1. Takes initial conditions as inputs.

            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the first stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx0MID = _mx0[spin], mx0LHS = _mx0[spinLHS], mx0RHS = _mx0[spinRHS];
            double my0MID = _my0[spin], my0LHS = _my0[spinLHS], my0RHS = _my0[spinRHS];
            double mz0MID = _mz0[spin], mz0LHS = _mz0[spinLHS], mz0RHS = _mz0[spinRHS];

            double hX0; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                if (_hasStaticDrive)
                    hX0 = _exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS + _dynamicBiasField;
                else if (!_hasStaticDrive)
                    hX0 = _exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS + _dynamicBiasField * cos(_drivingAngFreq * t0);
            } else
                // All spins along x which are not within the driving region
                hX0 = _exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS;

            double hY0 = _exchangeVec[spinLHS] * my0LHS + _exchangeVec[spin] * my0RHS;
            double hZ0 = _exchangeVec[spinLHS] * mz0LHS + _exchangeVec[spin] * mz0RHS + GV.GetStaticBiasField();

            double mxK1, myK1, mzK1; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK2.
                mxK1 = _gyroMagConst * (- (_gilbertVector[spin] * hY0 * mx0MID * my0MID) + hY0 * mz0MID - hZ0 * (my0MID + _gilbertVector[spin] * mx0MID * mz0MID) + _gilbertVector[spin] * hX0 * (pow(my0MID,2) + pow(mz0MID,2)));
                myK1 = _gyroMagConst * (-(hX0 * mz0MID) + hZ0 * (mx0MID - _gilbertVector[spin] * my0MID * mz0MID) + _gilbertVector[spin] * (hY0 * pow(mx0MID,2) - hX0 * mx0MID * my0MID + hY0 * pow(mz0MID,2)));
                mzK1 = _gyroMagConst * (hX0 * my0MID + _gilbertVector[spin] * hZ0 * (pow(mx0MID,2) + pow(my0MID,2)) - _gilbertVector[spin]*hX0*mx0MID*mz0MID - hY0 * (mx0MID + _gilbertVector[spin] * my0MID * mz0MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
                mxK1 = -1.0 * _gyroMagConst * (my0MID * hZ0 - mz0MID * hY0);
                myK1 =        _gyroMagConst * (mx0MID * hZ0 - mz0MID * hX0);
                mzK1 = -1.0 * _gyroMagConst * (mx0MID * hY0 - my0MID * hX0);
            }

            // Find (m0 + k1/2) for each spin, which is used in the next stage.
            mx1[spin] = mx0MID + mxK1 * _stepsizeHalf;
            my1[spin] = my0MID + myK1 * _stepsizeHalf;
            mz1[spin] = mz0MID + mzK1 * _stepsizeHalf;
        }
        // The estimations of the m-components' values for the next iteration.
        std::vector<double> mx2(GV.GetNumSpins() + 2,0), my2(GV.GetNumSpins() + 2,0), mz2(GV.GetNumSpins() + 2,0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            /* RK2 Step 2. Uses the previously found m1 values, as well as the initial conditions, to determine the
             * value of the slope (k2) at the midpoint. In RK2, the values of k1 and k2 can then be jointly used to
             * estimate the next point of the function through a weighted average of k1 & k2.
             */
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the second stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx1MID = mx1[spin], mx1LHS = mx1[spinLHS], mx1RHS = mx1[spinRHS];
            double my1MID = my1[spin], my1LHS = my1[spinLHS], my1RHS = my1[spinRHS];
            double mz1MID = mz1[spin], mz1LHS = mz1[spinLHS], mz1RHS = mz1[spinRHS];

            double hX1; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // If a spin is driven during Stage 1 of an RK method, then it must be driven throughout the rest of the method's stages. Note the different time value used
                if (_hasStaticDrive)
                    hX1 = _exchangeVec[spinLHS] * mx1LHS + _exchangeVec[spin] * mx1RHS + _dynamicBiasField;
                else if (!_hasStaticDrive)
                    hX1 = _exchangeVec[spinLHS] * mx1LHS + _exchangeVec[spin] * mx1RHS + _dynamicBiasField * cos(_drivingAngFreq * t0HalfStep);
            } else
                hX1 = _exchangeVec[spinLHS] * mx1LHS + _exchangeVec[spin] * mx1RHS;

            double hY1 = _exchangeVec[spinLHS] * my1LHS + _exchangeVec[spin] * my1RHS;
            double hZ1 = _exchangeVec[spinLHS] * mz1LHS + _exchangeVec[spin] * mz1RHS + GV.GetStaticBiasField();

            double mxK2, myK2, mzK2;
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation)
                mxK2 = _gyroMagConst * (- (_gilbertVector[spin] * hY1 * mx1MID * my1MID) + hY1 * mz1MID - hZ1 * (my1MID + _gilbertVector[spin] * mx1MID * mz1MID) + _gilbertVector[spin] * hX1 * (pow(my1MID,2) + pow(mz1MID,2)));
                myK2 = _gyroMagConst * (-(hX1 * mz1MID) + hZ1 * (mx1MID - _gilbertVector[spin] * my1MID * mz1MID) + _gilbertVector[spin] * (hY1 * pow(mx1MID,2) - hX1 * mx1MID * my1MID + hY1 * pow(mz1MID,2)));
                mzK2 = _gyroMagConst * (hX1 * my1MID + _gilbertVector[spin] * hZ1 * (pow(mx1MID,2) + pow(my1MID,2)) - _gilbertVector[spin]*hX1*mx1MID*mz1MID - hY1 * (mx1MID + _gilbertVector[spin] * my1MID * mz1MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation)
                mxK2 = -1.0 * _gyroMagConst * (my1MID * hZ1 - mz1MID * hY1);
                myK2 =        _gyroMagConst * (mx1MID * hZ1 - mz1MID * hX1);
                mzK2 = -1.0 * _gyroMagConst * (mx1MID * hY1 - my1MID * hX1);
            }

            mx2[spin] = _mx0[spin] + mxK2 * _stepsize;
            my2[spin] = _my0[spin] + myK2 * _stepsize;
            mz2[spin] = _mz0[spin] + mzK2 * _stepsize;

            if (_shouldTrackMValues) {
                double mIterationNorm = sqrt(pow(mx2[spin], 2) + pow(my2[spin], 2) + pow(mz2[spin], 2));
                if ((_largestMNorm) > (1.0 - mIterationNorm)) { _largestMNorm = (1.0 - mIterationNorm); }
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        _mx0.clear();
        _my0.clear();
        _mz0.clear();
        mx1.clear();
        my1.clear();
        mz1.clear();

        SaveDataToFile(mxRK2File, mx2, iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        _mx0 = mx2;
        _my0 = my2;
        _mz0 = mz2;

        if (iteration == _forceStopAtIteration)
            exit(0);

        _totalTime += _stepsize;
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

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

void Numerical_Methods_Class::RK2MidpointAFM() {

    progressbar bar(100);

    InformUserOfCodeType("RK2 Midpoint (AFM)");

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    //std::ofstream myRK2File(GV.GetFilePath() + "rk2_my_" + GV.GetFileNameBase() + ".csv");
    //std::ofstream mzRK2File(GV.GetFilePath() + "rk2_mz_" + GV.GetFileNameBase() + ".csv");
    CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
    //CreateFileHeader(myRK2File, "RK2 Midpoint (AFM)");
    //CreateFileHeader(mzRK2File, "RK2 Midpoint (AFM)");

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {

        if (_iterationEnd >= 100 && iteration % (_iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        TestShockwaveConditions(iteration);

        double t0 = _totalTime, t0HalfStep = _totalTime + _stepsizeHalf;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(GV.GetNumSpins() + 2, 0), my1(GV.GetNumSpins() + 2, 0), mz1(GV.GetNumSpins() + 2, 0);

        // Excludes the 0th and last spins as they will always be zero-valued (end, pinned spins)
        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            // RK2 Stage 1. Takes initial conditions as inputs.

            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the first stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx0MID = _mx0[spin], mx0LHS = _mx0[spinLHS], mx0RHS = _mx0[spinRHS];
            double my0MID = _my0[spin], my0LHS = _my0[spinLHS], my0RHS = _my0[spinRHS];
            double mz0MID = _mz0[spin], mz0LHS = _mz0[spinLHS], mz0RHS = _mz0[spinRHS];

            double hX0; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                if (_hasStaticDrive)
                    hX0 = -1.0 * (_exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS + _dynamicBiasField);
                else if (!_hasStaticDrive)
                    hX0 = -1.0 * (_exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS) + _dynamicBiasField * cos(_drivingAngFreq * t0);
            } else
                // All spins along x which are not within the driving region
                hX0 = -1.0 * (_exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS);

            double hY0 = -1.0 * (_exchangeVec[spinLHS] * my0LHS + _exchangeVec[spin] * my0RHS);

            //double hZ0 = GV.GetStaticBiasField() + _anisotropyField * mz0MID + (_exchangeVec[spinLHS] * mz0LHS + _exchangeVec[spin] * mz0RHS);
            double hZ0 = 0;
            if (mz0MID > 0)
                hZ0 = GV.GetStaticBiasField() + _anisotropyField - (_exchangeVec[spinLHS] * mz0LHS + _exchangeVec[spin] * mz0RHS);
            else if (mz0MID < 0)
                hZ0 = GV.GetStaticBiasField() - _anisotropyField - (_exchangeVec[spinLHS] * mz0LHS + _exchangeVec[spin] * mz0RHS);

            double mxK1, myK1, mzK1; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK2.
                mxK1 = _gyroMagConst * (- (_gilbertVector[spin] * hY0 * mx0MID * my0MID) + hY0 * mz0MID - hZ0 * (my0MID + _gilbertVector[spin] * mx0MID * mz0MID) + _gilbertVector[spin] * hX0 * (pow(my0MID,2) + pow(mz0MID,2)));
                myK1 = _gyroMagConst * (-(hX0 * mz0MID) + hZ0 * (mx0MID - _gilbertVector[spin] * my0MID * mz0MID) + _gilbertVector[spin] * (hY0 * pow(mx0MID,2) - hX0 * mx0MID * my0MID + hY0 * pow(mz0MID,2)));
                mzK1 = _gyroMagConst * (hX0 * my0MID + _gilbertVector[spin] * hZ0 * (pow(mx0MID,2) + pow(my0MID,2)) - _gilbertVector[spin]*hX0*mx0MID*mz0MID - hY0 * (mx0MID + _gilbertVector[spin] * my0MID * mz0MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
                mxK1 = -1.0 * _gyroMagConst * (my0MID * hZ0 - mz0MID * hY0);
                myK1 =        _gyroMagConst * (mx0MID * hZ0 - mz0MID * hX0);
                mzK1 = -1.0 * _gyroMagConst * (mx0MID * hY0 - my0MID * hX0);
            }

            mx1[spin] = mx0MID + mxK1 * _stepsizeHalf;
            my1[spin] = my0MID + myK1 * _stepsizeHalf;
            mz1[spin] = mz0MID + mzK1 * _stepsizeHalf;
        }
        // The estimations of the m-components' values for the next iteration.
        std::vector<double> mx2(GV.GetNumSpins() + 2,0), my2(GV.GetNumSpins() + 2,0), mz2(GV.GetNumSpins() + 2,0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            /*
             * RK2 Step 2. Uses the previously found m1 values, as well as the initial conditions, to determine the
             * value of the slope (k2) at the midpoint. In RK2, the values of k1 and k2 can then be jointly used to
             * estimate the next point of the function through a weighted average of k1 & k2.
             */
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the second stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx1MID = mx1[spin], mx2LHS = mx1[spinLHS], mx2RHS = mx1[spinRHS];
            double my1MID = my1[spin], my2LHS = my1[spinLHS], my2RHS = my1[spinRHS];
            double mz1MID = mz1[spin], mz2LHS = mz1[spinLHS], mz2RHS = mz1[spinRHS];

            // The effective field (H_eff) components acting upon each spin
            double hX1;
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // If a spin is driven during Stage 1 of an RK method, then it must be driven throughout the rest of the method's stages. Note the different time value used
                if (_hasStaticDrive)
                    hX1 = -1.0 * (_exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS + _dynamicBiasField);
                else if (!_hasStaticDrive)
                    hX1 = -1.0 * (_exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS) + _dynamicBiasField * cos(_drivingAngFreq * t0HalfStep);
            } else
                hX1 = -1.0 * (_exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS);

            double hY1 = -1.0 * (_exchangeVec[spinLHS] * my2LHS + _exchangeVec[spin] * my2RHS);

            double hZ1 = 0;
            if (mz1MID > 0)
                hZ1= GV.GetStaticBiasField() + _anisotropyField - (_exchangeVec[spinLHS] * mz2LHS + _exchangeVec[spin] * mz2RHS);
            else if (mz1MID < 0)
                hZ1= GV.GetStaticBiasField() - _anisotropyField - (_exchangeVec[spinLHS] * mz2LHS + _exchangeVec[spin] * mz2RHS);

            double mxK2, myK2, mzK2;
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation)
                mxK2 = _gyroMagConst * (- (_gilbertVector[spin] * hY1 * mx1MID * my1MID) + hY1 * mz1MID - hZ1 * (my1MID + _gilbertVector[spin] * mx1MID * mz1MID) + _gilbertVector[spin] * hX1 * (pow(my1MID,2) + pow(mz1MID,2)));
                myK2 = _gyroMagConst * (-(hX1 * mz1MID) + hZ1 * (mx1MID - _gilbertVector[spin] * my1MID * mz1MID) + _gilbertVector[spin] * (hY1 * pow(mx1MID,2) - hX1 * mx1MID * my1MID + hY1 * pow(mz1MID,2)));
                mzK2 = _gyroMagConst * (hX1 * my1MID + _gilbertVector[spin] * hZ1 * (pow(mx1MID,2) + pow(my1MID,2)) - _gilbertVector[spin]*hX1*mx1MID*mz1MID - hY1 * (mx1MID + _gilbertVector[spin] * my1MID * mz1MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation)
                mxK2 = -1.0 * _gyroMagConst * (my1MID * hZ1 - mz1MID * hY1);
                myK2 =        _gyroMagConst * (mx1MID * hZ1 - mz1MID * hX1);
                mzK2 = -1.0 * _gyroMagConst * (mx1MID * hY1 - my1MID * hX1);
            }

            mx2[spin] = _mx0[spin] + mxK2 * _stepsize;
            my2[spin] = _my0[spin] + myK2 * _stepsize;
            mz2[spin] = _mz0[spin] + mzK2 * _stepsize;

            if (_shouldTrackMValues) {
                double normOfIteration = sqrt(pow(mx2[spin], 2) + pow(my2[spin], 2) + pow(mz2[spin], 2));
                if ((_largestMNorm) > (1.0 - normOfIteration)) { _largestMNorm = (1.0 - normOfIteration); }
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        _mx0.clear();
        _my0.clear();
        _mz0.clear();
        mx1.clear();
        my1.clear();
        mz1.clear();

        SaveDataToFile(mxRK2File, mx2, iteration);
        //SaveDataToFile(myRK2File, my2, iteration);
        //SaveDataToFile(mzRK2File, mz2, iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        if (std::isnan(mx2[10]) or std::isnan(my2[10]) or std::isnan(mz2[10])) {
            std::cout << "\nNaN detected at iteration: " << iteration << std::endl; exit(1);
        }

        _mx0 = mx2;
        _my0 = my2;
        _mz0 = mz2;

        if (iteration == _forceStopAtIteration)
            exit(0);

        _totalTime += _stepsize;
    }
    // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();
    //myRK2File.close();
    //mzRK2File.close();

    if (_shouldTrackMValues)
        std::cout << "\nMax norm. value of M is: " << _largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::RK2MidpointFMForTesting() {

    progressbar bar(100);

    InformUserOfCodeType("RK2 Midpoint (Testing)");

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream SystemAllValuesPart1(GV.GetFilePath() + "rk2_HalfStep_" + GV.GetFileNameBase() + ".csv");
    std::ofstream SystemAllValuesPart2(GV.GetFilePath() + "rk2_FullStep_" + GV.GetFileNameBase() + ".csv");

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {

        if (iteration == _forceStopAtIteration) {
            SystemAllValuesPart1 << "Spin,Iteration,t0,mx0,my0,mz0,hX0,hY0,hZ0,mxK1,myK1,mzK1,mx1Err,my1Err,mz1Err,mx1,my1,mz1\n";
            SystemAllValuesPart2 << "Spin,Iteration,t0h,mx1,my1,mz1,hX1,hY1,hZ1,mxK2,myK2,mzK2,mx2Err,my2Err,mz2Err,mx2,my2,mz2,M Norm, Delta Norm\n";
        }
        if (_iterationEnd >= 100 && iteration % (_iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        TestShockwaveConditions(iteration);
        _totalTime += _stepsize;
        double t0 = _totalTime, t0HalfStep = _totalTime + _stepsizeHalf;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(GV.GetNumSpins() + 2, 0), my1(GV.GetNumSpins() + 2, 0), mz1(GV.GetNumSpins() + 2, 0);

        // Excludes the 0th and last spins as they will always be zero-valued (end, pinned spins)
        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            /*
             * The first stage is based upon finding the value of the slope at the beginning of the interval (k1). This
             * stage takes the start conditions as an input.
             */
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the first stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx0MID = _mx0[spin], mx0LHS = _mx0[spinLHS], mx0RHS = _mx0[spinRHS];
            double my0MID = _my0[spin], my0LHS = _my0[spinLHS], my0RHS = _my0[spinRHS];
            double mz0MID = _mz0[spin], mz0LHS = _mz0[spinLHS], mz0RHS = _mz0[spinRHS];

            double hX0; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                if (_hasStaticDrive)
                    hX0 = _exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS + _dynamicBiasField;
                else if (!_hasStaticDrive)
                    hX0 = _exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS + _dynamicBiasField * cos(_drivingAngFreq * t0);
            } else
                // All spins along x which are not within the driving region
                hX0 = _exchangeVec[spinLHS] * mx0LHS + _exchangeVec[spin] * mx0RHS;

            // No changes are made to the effective field in the y-direction
            double hY0 = _exchangeVec[spinLHS] * my0LHS + _exchangeVec[spin] * my0RHS;
            // The static bias field is applied in the z-direction
            double hZ0 = _exchangeVec[spinLHS] * mz0LHS + _exchangeVec[spin] * mz0RHS + GV.GetStaticBiasField();

            double mxK1, myK1, mzK1; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK2.

                mxK1 = _gyroMagConst * (- (_gilbertVector[spin] * hY0 * mx0MID * my0MID) + hY0 * mz0MID - hZ0 * (my0MID + _gilbertVector[spin]*mx0MID*mz0MID) + _gilbertVector[spin] * hX0 * (pow(my0MID,2) + pow(mz0MID,2)));
                myK1 = _gyroMagConst * (-(hX0*mz0MID) + hZ0 * (mx0MID - _gilbertVector[spin] * my0MID * mz0MID) + _gilbertVector[spin] * (hY0 * pow(mx0MID,2) - hX0 * mx0MID * my0MID + hY0 * pow(mz0MID,2)));
                mzK1 = _gyroMagConst * (hX0 * my0MID + _gilbertVector[spin] * hZ0*(pow(mx0MID,2) + pow(my0MID,2)) - _gilbertVector[spin]*hX0*mx0MID*mz0MID - hY0 * (mx0MID + _gilbertVector[spin] * my0MID * mz0MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
                mxK1 = -1 * _gyroMagConst * (my0MID * hZ0 - mz0MID * hY0);
                myK1 = _gyroMagConst * (mx0MID * hZ0 - mz0MID * hX0);
                mzK1 = -1 * _gyroMagConst * (mx0MID * hY0 - my0MID * hX0);
            }

            // Find (m0 + k1/2) for each spin, which is used in the next stage.
            mx1[spin] = mx0MID + mxK1 * _stepsizeHalf;
            my1[spin] = my0MID + myK1 * _stepsizeHalf;
            mz1[spin] = mz0MID + mzK1 * _stepsizeHalf;

            double mx1err, my1err, mz1err;
            mx1err = -1.0 * ((hY0 * mx0MID *  my0MID - hX0 * pow(my0MID,2) + hZ0 * mx0MID * mz0MID - hX0 * pow(mz0MID, 2)) * _gilbertVector[spin]);
            my1err = -1.0 * ((-hY0 * pow(mx0MID,2) + hX0 * mx0MID * my0MID + hZ0 * my0MID * mz0MID - hY0 * pow(mz0MID,2)) * _gilbertVector[spin]);
            mz1err = -1.0 * ((-hZ0 * pow(mx0MID, 2) - hZ0 * pow(my0MID, 2) + hX0 * mx0MID * mz0MID + hY0 * my0MID * mz0MID) * _gilbertVector[spin]);

            if (iteration == _forceStopAtIteration)
                SystemAllValuesPart1 << spin << ", " << iteration + 1 << ", " << _totalTime << ", " << mx0MID << ", " << my0MID << ", " << mz0MID << ", " << hX0 << ", " << hY0 << ", " << hZ0 << ", " << mxK1 << ", " << myK1 << ", " << mzK1 << ", " << mx1err << ", " << my1err << ", " << mz1err << ", " << mx1[spin] << ", " << my1[spin] << ", " << mz1[spin] << "\n";
        }
        if (iteration == _forceStopAtIteration)
            SystemAllValuesPart1 << "\n\n\n\n\n\n\n\n\n\n";

        // The estimations of the m-components' values for the next iteration.
        std::vector<double> mx2(GV.GetNumSpins() + 2,0), my2(GV.GetNumSpins() + 2,0), mz2(GV.GetNumSpins() + 2,0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            /** 
             * The second stage uses the previously found m1 values, as well as the initial conditions, to determine the
             * value of the slope (k2) at the midpoint. In RK2, the values of k1 and k2 can then be jointly used to
             * estimate the next point of the function through a weighted average of k1 & k2.
             */
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the second stage for the: current spin site (CUR); site to the left (LHS); site to the right (RHS)
            double mx1MID = mx1[spin], mx2LHS = mx1[spinLHS], mx2RHS = mx1[spinRHS];
            double my1MID = my1[spin], my2LHS = my1[spinLHS], my2RHS = my1[spinRHS];
            double mz1MID = mz1[spin], mz2LHS = mz1[spinLHS], mz2RHS = mz1[spinRHS];

            double hX1; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // If a spin is driven during Stage 1 of an RK method, then it must be driven throughout the rest of the method's stages. Note the different time value used
                if (_hasStaticDrive)
                    hX1 = _exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS + _dynamicBiasField;
                else if (!_hasStaticDrive)
                    hX1 = _exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS + _dynamicBiasField * cos(_drivingAngFreq * t0HalfStep);
            } else
                hX1 = _exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS;

            double hY1 = _exchangeVec[spinLHS] * my2LHS + _exchangeVec[spin] * my2RHS;
            double hZ1 = _exchangeVec[spinLHS] * mz2LHS + _exchangeVec[spin] * mz2RHS + GV.GetStaticBiasField();

            double mxK2, myK2, mzK2;
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation)
                mxK2 = _gyroMagConst * (- (_gilbertVector[spin] * hY1 * mx1MID * my1MID) + hY1 * mz1MID - hZ1 * (my1MID + _gilbertVector[spin]*mx1MID*mz1MID) + _gilbertVector[spin] * hX1 * (pow(my1MID,2) + pow(mz1MID,2)));
                myK2 = _gyroMagConst * (-(hX1*mz1MID) + hZ1 * (mx1MID - _gilbertVector[spin] * my1MID * mz1MID) + _gilbertVector[spin] * (hY1 * pow(mx1MID,2) - hX1 * mx1MID * my1MID + hY1 * pow(mz1MID,2)));
                mzK2 = _gyroMagConst * (hX1 * my1MID + _gilbertVector[spin] * hZ1*(pow(mx1MID,2) + pow(my1MID,2)) - _gilbertVector[spin]*hX1*mx1MID*mz1MID - hY1 * (mx1MID + _gilbertVector[spin] * my1MID * mz1MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation)
                mxK2 = -1 * _gyroMagConst * (my1MID * hZ1 - mz1MID * hY1);
                myK2 =      _gyroMagConst * (mx1MID * hZ1 - mz1MID * hX1);
                mzK2 = -1 * _gyroMagConst * (mx1MID * hY1 - my1MID * hX1);
            }

            double mx2err, my2err, mz2err;
            mx2err = -1.0 * ((hY1 * mx1MID *  my1MID - hX1 * pow(my1MID,2) + hZ1 * mx1MID * mz1MID - hX1 * pow(mz1MID, 2)) * _gilbertVector[spin]);
            my2err = -1.0 * ((-hY1 * pow(mx1MID,2) + hX1 * mx1MID * my1MID + hZ1 * my1MID * mz1MID - hY1 * pow(mz1MID,2)) * _gilbertVector[spin]);
            mz2err = -1.0 * ((-hZ1 * pow(mx1MID, 2) - hZ1 * pow(my1MID, 2) + hX1 * mx1MID * mz1MID + hY1 * my1MID * mz1MID) * _gilbertVector[spin]);

            mx2[spin] = _mx0[spin] + mxK2 * _stepsize;
            my2[spin] = _my0[spin] + myK2 * _stepsize;
            mz2[spin] = _mz0[spin] + mzK2 * _stepsize;

            double mSumTotal = mx2[spin] + my2[spin] + mz2[spin];
            double mSumNorm = sqrt(pow(mx2[spin], 2) + pow(my2[spin], 2) + pow(mz2[spin], 2));

            if (iteration == _forceStopAtIteration)
                SystemAllValuesPart2 << spin << ", " << iteration + 1 << ", " << t0HalfStep << ", " << mx1MID << ", " << my1MID << ", " << mz1MID << ", " << hX1 << ", " << hY1 << ", " << hZ1 << ", " << mxK2 << ", " << myK2 << ", " << mzK2 << ", " << mx2err << "," << my2err << ", " << mz2err << ", " <<mx2[spin] << ", " << my2[spin] << ", " << mz2[spin] << ", " << mSumNorm << ", " << mSumNorm - 1.0 << "\n";

        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.
        if (iteration == _forceStopAtIteration)
            SystemAllValuesPart2 << "\n\n\n\n\n\n\n\n\n\n";
        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        _mx0.clear();
        _my0.clear();
        _mz0.clear();
        mx1.clear();
        my1.clear();
        mz1.clear();

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        _mx0 = mx2;
        _my0 = my2;
        _mz0 = mz2;

        if (iteration == _forceStopAtIteration) {
            std::cout << "Force Stop activated at iteration #" << iteration;
            exit(0);
        }

        _totalTime += _stepsize; // Bob has stepsize increment at the start of the loop
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    SystemAllValuesPart1.close();
    SystemAllValuesPart2.close();

    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::RK4MidpointFM() {

    progressbar bar(100);

    InformUserOfCodeType("RK4 Midpoint");

    //########################################################################################################################

    std::ofstream mxRK4File(GV.GetFilePath() + "rk4_mx_" + GV.GetFileNameBase() + ".csv");
    CreateFileHeader(mxRK4File, "RK4 Midpoint");

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {

        double t0 = _totalTime, t0Half = _totalTime + _stepsizeHalf, t0h = _totalTime + _stepsize;

        if (_iterationEnd >= 100 && iteration % (_iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        TestShockwaveConditions(iteration);

        std::vector<double> mx1(GV.GetNumSpins()+2, 0), my1(GV.GetNumSpins()+2, 0), mz1(GV.GetNumSpins()+2, 0);
        std::vector<double> mxK1Vec (GV.GetNumSpins() + 2, 0), myK1Vec (GV.GetNumSpins() + 2, 0), mzK1Vec (GV.GetNumSpins() + 2, 0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            // RK4 Step 1
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the first stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx1MID = _mx0[spin], mx1LHS = _mx0[spinLHS], mx1RHS = _mx0[spinRHS];
            double my1MID = _my0[spin], my1LHS = _my0[spinLHS], my1RHS = _my0[spinRHS];
            double mz1MID = _mz0[spin], mz1LHS = _mz0[spinLHS], mz1RHS = _mz0[spinRHS];

            double hX1; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS)  {
                if (_hasStaticDrive)
                    hX1 = _exchangeVec[spinLHS] * mx1LHS + _exchangeVec[spin] * mx1RHS + _dynamicBiasField;
                else if (!_hasStaticDrive)
                    hX1 = _exchangeVec[spinLHS] * mx1LHS + _exchangeVec[spin] * mx1RHS + _dynamicBiasField * cos(_drivingAngFreq * t0);
            } else {
                hX1 = _exchangeVec[spinLHS] * mx1LHS + _exchangeVec[spin] * mx1RHS;
            }
            double hY1 = _exchangeVec[spinLHS] * my1LHS + _exchangeVec[spin] * my1RHS;
            double hZ1 = _exchangeVec[spinLHS] * mz1LHS + _exchangeVec[spin] * mz1RHS + GV.GetStaticBiasField();

            double mxK1, myK1, mzK1; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK4.
                mxK1 = _gyroMagConst * (- (_gilbertVector[spin] * hY1 * mx1MID * my1MID) + hY1 * mz1MID - hZ1 * (my1MID + _gilbertVector[spin] * mx1MID * mz1MID) + _gilbertVector[spin] * hX1 * (pow(my1MID,2) + pow(mz1MID,2)));
                myK1 = _gyroMagConst * (-(hX1 * mz1MID) + hZ1 * (mx1MID - _gilbertVector[spin] * my1MID * mz1MID) + _gilbertVector[spin] * (hY1 * pow(mx1MID,2) - hX1 * mx1MID * my1MID + hY1 * pow(mz1MID,2)));
                mzK1 = _gyroMagConst * (hX1 * my1MID + _gilbertVector[spin] * hZ1 * (pow(mx1MID,2) + pow(my1MID,2)) - _gilbertVector[spin] * hX1 * mx1MID * mz1MID - hY1 * (mx1MID + _gilbertVector[spin] * my1MID * mz1MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK4.
                mxK1 = -1.0 * _gyroMagConst * (my1MID * hZ1 - mz1MID * hY1);
                myK1 =        _gyroMagConst * (mx1MID * hZ1 - mz1MID * hX1);
                mzK1 = -1.0 * _gyroMagConst * (mx1MID * hY1 - my1MID * hX1);
            }

            mxK1Vec[spin] = mxK1;
            myK1Vec[spin] = myK1;
            mzK1Vec[spin] = mzK1;

            // Find (m0 + k1/2) for each spin, which is used in the next stage.
            mx1[spin] = mx1MID + mxK1 * _stepsizeHalf;
            my1[spin] = my1MID + myK1 * _stepsizeHalf;
            mz1[spin] = mz1MID + mzK1 * _stepsizeHalf;
        }

        std::vector<double> mx2 (GV.GetNumSpins() + 2, 0), my2 (GV.GetNumSpins() + 2, 0), mz2 (GV.GetNumSpins() + 2, 0);
        std::vector<double> mxK2Vec (GV.GetNumSpins() + 2, 0), myK2Vec (GV.GetNumSpins() + 2, 0), mzK2Vec (GV.GetNumSpins() + 2, 0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            // RK4 Step 2
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the second stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx2MID = mx1[spin], mx2LHS = mx1[spinLHS], mx2RHS = mx1[spinRHS];
            double my2MID = my1[spin], my2LHS = my1[spinLHS], my2RHS = my1[spinRHS];
            double mz2MID = mz1[spin], mz2LHS = mz1[spinLHS], mz2RHS = mz1[spinRHS];

            double hX2; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS)  {
                if (_hasStaticDrive)
                    hX2 = _exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS + _dynamicBiasField;
                else if (!_hasStaticDrive)
                    hX2 = _exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS + _dynamicBiasField * cos(_drivingAngFreq * t0Half);
            } else {
                hX2 = _exchangeVec[spinLHS] * mx2LHS + _exchangeVec[spin] * mx2RHS;
            }
            double hY2 = _exchangeVec[spinLHS] * my2LHS + _exchangeVec[spin] * my2RHS;
            double hZ2 = _exchangeVec[spinLHS] * mz2LHS + _exchangeVec[spin] * mz2RHS + GV.GetStaticBiasField();

            double mxK2, myK2, mzK2; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the second stage of RK4.
                mxK2 = _gyroMagConst * (- (_gilbertVector[spin] * hY2 * mx2MID * my2MID) + hY2 * mz2MID - hZ2 * (my2MID + _gilbertVector[spin] * mx2MID * mz2MID) + _gilbertVector[spin] * hX2 * (pow(my2MID,2) + pow(mz2MID,2)));
                myK2 = _gyroMagConst * (-(hX2 * mz2MID) + hZ2 * (mx2MID - _gilbertVector[spin] * my2MID * mz2MID) + _gilbertVector[spin] * (hY2 * pow(mx2MID,2) - hX2 * mx2MID * my2MID + hY2 * pow(mz2MID,2)));
                mzK2 = _gyroMagConst * (hX2 * my2MID + _gilbertVector[spin] * hZ2 * (pow(mx2MID,2) + pow(my2MID,2)) - _gilbertVector[spin] * hX2 * mx2MID * mz2MID - hY2 * (mx2MID + _gilbertVector[spin] * my2MID * mz2MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the second stage of RK4.
                mxK2 = -1.0 * _gyroMagConst * (my2MID * hZ2 - mz2MID * hY2);
                myK2 =        _gyroMagConst * (mx2MID * hZ2 - mz2MID * hX2);
                mzK2 = -1.0 * _gyroMagConst * (mx2MID * hY2 - my2MID * hX2);
            }

            mxK2Vec[spin] = mxK2;
            myK2Vec[spin] = myK2;
            mzK2Vec[spin] = mzK2;

            // Find (m0 + k2/2) for each spin, which is used in the next stage.
            mx2[spin] = _mx0[spin] + mxK2 * _stepsizeHalf;
            my2[spin] = _my0[spin] + myK2 * _stepsizeHalf;
            mz2[spin] = _mz0[spin] + mzK2 * _stepsizeHalf;
        }

        mx1.clear(); my1.clear(); mz1.clear(); // No longer required so memory can be freed
        std::vector<double> mx3 (GV.GetNumSpins() + 2, 0), my3 (GV.GetNumSpins() + 2, 0), mz3 (GV.GetNumSpins() + 2, 0);
        std::vector<double> mxK3Vec (GV.GetNumSpins() + 2, 0), myK3Vec (GV.GetNumSpins() + 2, 0), mzK3Vec (GV.GetNumSpins() + 2, 0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            // RK4 Step 3
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the third stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx3MID = mx2[spin], mx3LHS = mx2[spinLHS], mx3RHS = mx2[spinRHS];
            double my3MID = my2[spin], my3LHS = my2[spinLHS], my3RHS = my2[spinRHS];
            double mz3MID = mz2[spin], mz3LHS = mz2[spinLHS], mz3RHS = mz2[spinRHS];

            double hX3; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS)  {
                if (_hasStaticDrive)
                    hX3 = _exchangeVec[spinLHS] * mx3LHS + _exchangeVec[spin] * mx3RHS + _dynamicBiasField;
                else if (!_hasStaticDrive)
                    hX3 = _exchangeVec[spinLHS] * mx3LHS + _exchangeVec[spin] * mx3RHS + _dynamicBiasField * cos(_drivingAngFreq * t0Half);
            } else {
                hX3 = _exchangeVec[spinLHS] * mx3LHS + _exchangeVec[spin] * mx3RHS;
            }
            double hY3 = _exchangeVec[spinLHS] * my3LHS + _exchangeVec[spin] * my3RHS;
            double hZ3 = _exchangeVec[spinLHS] * mz3LHS + _exchangeVec[spin] * mz3RHS + GV.GetStaticBiasField();

            double mxK3, myK3, mzK3; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK4.
                mxK3 = _gyroMagConst * (- (_gilbertVector[spin] * hY3 * mx3MID * my3MID) + hY3 * mz3MID - hZ3 * (my3MID + _gilbertVector[spin] * mx3MID * mz3MID) + _gilbertVector[spin] * hX3 * (pow(my3MID,2) + pow(mz3MID,2)));
                myK3 = _gyroMagConst * (-(hX3 * mz3MID) + hZ3 * (mx3MID - _gilbertVector[spin] * my3MID * mz3MID) + _gilbertVector[spin] * (hY3 * pow(mx3MID,2) - hX3 * mx3MID * my3MID + hY3 * pow(mz3MID,2)));
                mzK3 = _gyroMagConst * (hX3 * my3MID + _gilbertVector[spin] * hZ3 * (pow(mx3MID,2) + pow(my3MID,2)) - _gilbertVector[spin] * hX3 * mx3MID * mz3MID - hY3 * (mx3MID + _gilbertVector[spin] * my3MID * mz3MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK4.
                mxK3 = -1.0 * _gyroMagConst * (my3MID * hZ3 - mz3MID * hY3);
                myK3 =        _gyroMagConst * (mx3MID * hZ3 - mz3MID * hX3);
                mzK3 = -1.0 * _gyroMagConst * (mx3MID * hY3 - my3MID * hX3);
            }

            mxK3Vec[spin] = mxK3;
            myK3Vec[spin] = myK3;
            mzK3Vec[spin] = mzK3;

            // Find (m0 + k3) for each spin, which is used in the next stage.
            mx3[spin] = _mx0[spin] + mxK3 * _stepsize;
            my3[spin] = _my0[spin] + myK3 * _stepsize;
            mz3[spin] = _mz0[spin] + mzK3 * _stepsize;
        }

        mx2.clear(); my2.clear(); mz2.clear(); // No longer required so memory can be freed
        std::vector<double> mx4 (GV.GetNumSpins() + 2, 0), my4 (GV.GetNumSpins() + 2, 0), mz4 (GV.GetNumSpins() + 2, 0);
        std::vector<double> mxK (GV.GetNumSpins() + 2, 0), myK (GV.GetNumSpins() + 2, 0), mzK (GV.GetNumSpins() + 2, 0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            // RK4 Step 4
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the fourth stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx4MID = mx3[spin], mx4LHS = mx3[spinLHS], mx4RHS = mx3[spinRHS];
            double my4MID = my3[spin], my4LHS = my3[spinLHS], my4RHS = my3[spinRHS];
            double mz4MID = mz3[spin], mz4LHS = mz3[spinLHS], mz4RHS = mz3[spinRHS];

            double hX4; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS)  {
                if (_hasStaticDrive)
                    hX4 = _exchangeVec[spinLHS] * mx4LHS + _exchangeVec[spin] * mx4RHS + _dynamicBiasField;
                if (!_hasStaticDrive)
                    hX4 = _exchangeVec[spinLHS] * mx4LHS + _exchangeVec[spin] * mx4RHS + _dynamicBiasField * cos(_drivingAngFreq * t0h);

            } else {
                hX4 = _exchangeVec[spinLHS] * mx4LHS + _exchangeVec[spin] * mx4RHS;
            }
            double hY4 = _exchangeVec[spinLHS] * my4LHS + _exchangeVec[spin] * my4RHS;
            double hZ4 = _exchangeVec[spinLHS] * mz4LHS + _exchangeVec[spin] * mz4RHS + GV.GetStaticBiasField();

            double mxK4, myK4, mzK4; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the first stage of RK4.
                mxK4 = _gyroMagConst * (- (_gilbertVector[spin] * hY4 * mx4MID * my4MID) + hY4 * mz4MID - hZ4 * (my4MID + _gilbertVector[spin] * mx4MID * mz4MID) + _gilbertVector[spin] * hX4 * (pow(my4MID,2) + pow(mz4MID,2)));
                myK4 = _gyroMagConst * (-(hX4 * mz4MID) + hZ4 * (mx4MID - _gilbertVector[spin] * my4MID * mz4MID) + _gilbertVector[spin] * (hY4 * pow(mx4MID,2) - hX4 * mx4MID * my4MID + hY4 * pow(mz4MID,2)));
                mzK4 = _gyroMagConst * (hX4 * my4MID + _gilbertVector[spin] * hZ4 * (pow(mx4MID,2) + pow(my4MID,2)) - _gilbertVector[spin] * hX4 * mx4MID * mz4MID - hY4 * (mx4MID + _gilbertVector[spin] * my4MID * mz4MID));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK4.
                mxK4 = -1.0 * _gyroMagConst * (my4MID * hZ4 - mz4MID * hY4);
                myK4 =        _gyroMagConst * (mx4MID * hZ4 - mz4MID * hX4);
                mzK4 = -1.0 * _gyroMagConst * (mx4MID * hY4 - my4MID * hX4);
            }

            mxK[spin] = (mxK1Vec[spin] + 2.0 * mxK2Vec[spin] + 2.0 * mxK3Vec[spin] + mxK4) * (_stepsize / 6.0);
            myK[spin] = (myK1Vec[spin] + 2.0 * myK2Vec[spin] + 2.0 * myK3Vec[spin] + myK4) * (_stepsize / 6.0);
            mzK[spin] = (mzK1Vec[spin] + 2.0 * mzK2Vec[spin] + 2.0 * mzK3Vec[spin] + mzK4) * (_stepsize / 6.0);

            // Find (m0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6) for each spin, as per the Runge-Kutta formula
            mx4[spin] = _mx0[spin] + mxK[spin];
            my4[spin] = _my0[spin] + myK[spin];
            mz4[spin] = _mz0[spin] + mzK[spin];

            if (_shouldTrackMValues) {
                double mIterationNorm = sqrt(pow(mx4[spin], 2) + pow(my4[spin], 2) + pow(mz4[spin], 2));
                if (_largestMNorm < (1.0 - mIterationNorm)) { _largestMNorm = mIterationNorm; }
            }

        }
        // Everything below here is part of the class method, but not the internal RK4 stage loops.

        // Clear all remaining vectors not required for writing processes
        mx3.clear(); my3.clear(); mz3.clear();
        mxK1Vec.clear(); myK1Vec.clear(); mzK1Vec.clear();
        mxK2Vec.clear(); myK2Vec.clear(); mzK2Vec.clear();
        mxK3Vec.clear(); myK3Vec.clear(); mzK3Vec.clear();

        SaveDataToFile(mxRK4File, mx4, iteration);
        // SaveDataToFile(_saveAllSpins, myRK4File, my4, iteration, _printFixedLines);
        // SaveDataToFile(_saveAllSpins, mzRK4File, mz4, iteration, _printFixedLines);

        // Set final value of current iteration as the starting value for the next iteration
        _mx0 = mx4;
        _my0 = my4;
        _mz0 = mz4;

        _totalTime += _stepsize;
    }
    // Final line of RK2 solver for all iterations. Everything below here occurs after RK4 method is complete.

    // Ensure files are closed in case of writing error.
    mxRK4File.close();
    // myRK4File.close();
    // mzRK4File.close();

    if (_shouldTrackMValues)
        std::cout << "\nMax norm. value of M is: " << _largestMNorm << std::endl;
    // Filename can be copy/pasted from C++ console to a Python console
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
    if (_saveAllSpins or _printFixedLines) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for (int i = 1; i <= GV.GetNumSpins(); i++) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if (_printFixedSites) {

        outputFileName << "Time" << ", ";
        for (int & fixed_out_val : _fixed_output_sites)
            std::cout << fixed_out_val << ", ";
        outputFileName << std::endl;

        //outputFileName << "Time" << ", "
        //               << static_cast<int>(14000) << ","
        //               << static_cast<int>(16000) << ","
        //               << static_cast<int>(18000) << ","
        //               << static_cast<int>(20000) << std::endl;

    } else {
        outputFileName << "Time" << ", "
                       << _drivingRegionLHS << ","
                       << static_cast<int>(_drivingRegionWidth / 2.0) << ","
                       << _drivingRegionRHS << ","
                       << static_cast<int>(GV.GetNumSpins() / 4.0) << ","
                       << static_cast<int>(GV.GetNumSpins() / 2.0) << ","
                       << static_cast<int>(3.0 * GV.GetNumSpins() / 4.0) << ","
                       << GV.GetNumSpins() << std::endl;
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
            outputFileName << (iteration * _stepsize) << ", ";
            for (int & fixed_out_val : _fixed_output_sites)
                std::cout << arrayToWrite[fixed_out_val] << ", ";
            outputFileName << std::endl;
            
            return;
        }
    }

    if (_saveAllSpins) {
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
        if (_saveAllSpins) {
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