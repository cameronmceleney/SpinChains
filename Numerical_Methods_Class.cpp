#include "Numerical_Methods_Class.h"

void Numerical_Methods_Class::NMSetup() {

    // ###################### Flags ######################
    _hasShockwave = false;
    _lhsDrive = true;
    _useLLG = true;
    _shouldTrackMValues = false;

    // ###################### Core Parameters ######################
    _biasFieldDriving = 3e-3;
    _drivingFreq = 2.5 * 1e9;
    _stepsize = 4e-15;
    _iterationEnd = static_cast<int>(1e3);

    // ###################### Shockwave Parameters ######################
    _iterToBeginShockwave = 0.0;
    _shockwaveScaling = 0;
    _shockwaveInit = _biasFieldDriving;
    _shockwaveMax = _shockwaveInit * _shockwaveScaling;
    _shockwaveGradientTime = 50e3;
    _shockwaveStepsize = (_shockwaveMax - _shockwaveInit) / _shockwaveGradientTime;

    // ###################### Data Output Parameters ######################
    _onlyShowFinalState = true;
    _saveAllSpins = false;
    _fixedPoints = false;
    _numberOfDataPoints = 100;

    // ###################### Damping Factors ######################
    _gilbertLower = _gilbertConst;
    _gilbertUpper = 1.0;

    // ###################### SpinChain Length Parameters ######################
    _drivingRegionWidth = 10;
    _numSpinsDamped = 0;
    _numSpinsInChain = GV.GetNumSpins();

    // ###################### Computations based upon other inputs ######################
    _drivingAngFreq = 2 * M_PI * _drivingFreq;
    _maxSimTime = _stepsize * _iterationEnd;
    _numberOfSpinPairs = GV.GetNumSpins() - 1;
    _stepsizeHalf = _stepsize / 2.0;
    GV.SetNumSpins(_numSpinsInChain + 2 * _numSpinsDamped);

    // ###################### Method Invocations ######################
    SetDrivingRegion(_lhsDrive); // Create driving region
    SetupVectors(); // Generate exchange and damping-region vectors needed for all RK methods.
}
void Numerical_Methods_Class::SetDrivingRegion(bool &useLHSDrive) {
    /**
     * Set up driving regions for the system. The LHS option is solely for drives from the left of the system. The RHS options contains the
     * drive from the right, as well as an option to drive from the centre.
     */
    if (useLHSDrive) {
        //Drives from the LHS, starting at _drivingRegionLHS
        _drivingRegionLHS = _numSpinsDamped + 1;  // The +1 is to correct the offset of adding a zeroth spin
        // _drivingRegionWidth = static_cast<int>(_numSpinsInChain * _regionScaling);
        _drivingRegionRHS = _drivingRegionLHS + _drivingRegionWidth;
    }
    else {
        // Drives from the RHS, starting at _drivingRegionRHS
        _drivingRegionRHS = GV.GetNumSpins() - _numSpinsDamped;
        // _drivingRegionWidth = static_cast<int>(_numSpinsInChain * _regionScaling);
        //_drivingRegionRHS = (_numSpinsInChain/2) +_numSpinsDamped + (_drivingRegionWidth / 2); // use for central drive
        _drivingRegionLHS = _drivingRegionRHS - _drivingRegionWidth + 1;  // The +1 is to correct the offset of adding a zeroth spin
    }
}
void Numerical_Methods_Class::SetupVectors() {

    SetupVectorsExchange();
    // SetupVectorsGilbert();
    GilbertVectorsBothSides();
}
void Numerical_Methods_Class::SetupVectorsExchange() {
    /*
     * Create the arrays which house the exchange integral values. There are options to have a non-uniform exchange coded in, as well as the option to
     * induce a 'kick' into the system by initialising certain spins to have differing parameters to their neighbours.
     */
    LinspaceClass SpinChainExchange;

    // The linearly spaced vector is saved as the class member '_chainJVals' simply to increase code readability
    SpinChainExchange.set_values(GV.GetExchangeMinVal(), GV.GetExchangeMaxVal(), _numberOfSpinPairs, true, true);
    _chainJVals = SpinChainExchange.generate_array();

    /*
    for (double i: _chainJVals) {
        std::cout << i << ", ";
    }
    exit(0);
     */

    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mxInitCond(GV.GetNumSpins(), _mxInit), myInitCond(GV.GetNumSpins(), _myInit), mzInitCond(GV.GetNumSpins(), _mzInit);
    // mxInitCond[0] = _mxInit; // Only perturb initial spin

    //int spinsEitherSide = 0;
    /*
    for (int i=0; i < 10; i++)
    {
        //mxInitCond[i] = 0.0001; // 0.00001 and 0.99999
        //myInitCond[i] = 0.0;
        //mzInitCond[i] = 0.99999;

        mxInitCond[i] = -0.0001;
        myInitCond[i] = 0.0;
        mzInitCond[i] = 0.99999;
    }*/


    // Appends initial conditions to the vectors
    _mx0.insert(_mx0.end(), mxInitCond.begin(), mxInitCond.end());
    _my0.insert(_my0.end(), myInitCond.begin(), myInitCond.end());
    _mz0.insert(_mz0.end(), mzInitCond.begin(), mzInitCond.end());

    // This zero is the (N+1)th spin on the RHS of the chain
    _mx0.push_back(0);
    _my0.push_back(0);
    _mz0.push_back(0);
}
void Numerical_Methods_Class::GilbertVectorsBothSides() {
    /*
     * Generate the damping regions that are appended to either end of the spin chain.
     */
    LinspaceClass GilbertDampingLHS;
    LinspaceClass GilbertDampingRHS;

    if (_numSpinsDamped < 0) {
        // Guard clause.
        std::cout << "numGilbert is less than zero!";
        exit(0);
    }

    std::vector<double> gilbertChain(_numSpinsInChain, _gilbertConst);

    GilbertDampingLHS.set_values(_gilbertUpper, _gilbertLower, _numSpinsDamped, true, false);
    GilbertDampingRHS.set_values(_gilbertLower, _gilbertUpper, _numSpinsDamped, true, false);
    std::vector<double> tempGilbertLHS = GilbertDampingLHS.generate_array();
    std::vector<double> tempGilbertRHS = GilbertDampingRHS.generate_array();

    _gilbertVector.insert(_gilbertVector.end(), tempGilbertLHS.begin(), tempGilbertLHS.end());
    _gilbertVector.insert(_gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());
    _gilbertVector.insert(_gilbertVector.end(), tempGilbertRHS.begin(), tempGilbertRHS.end());
    _gilbertVector.push_back(0);

    /*
    int count = 0;
    for (double i: _gilbertVector) {
        if (++count % 10 == 0)
            std::cout << std::setw(12) << i << std::endl;
        else
            std::cout << std::setw(12) << i << ", ";
    }
    exit(0);
     */

}

void Numerical_Methods_Class::RK2Original() {
    
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
                HeffX1K1 = _chainJVals[LHS_spin] * mx1LHS + _chainJVals[spin] * mx1RHS + _biasFieldDriving*cos(_drivingAngFreq * t0);
            } else {
                // The else statement includes all spins along x which are not within the driving region
                HeffX1K1 = _chainJVals[LHS_spin] * mx1LHS + _chainJVals[spin] * mx1RHS;
            }

            // No changes are made to the effective field in the y-direction
            HeffY1K1 = _chainJVals[LHS_spin] * my1LHS + _chainJVals[spin] * my1RHS;
            // The bias field is applied in the z-direction and so it contributes to the effective field in the z-direction
            HeffZ1K1 = _chainJVals[LHS_spin] * mz1LHS + _chainJVals[spin] * mz1RHS + GV.GetBiasField();

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
                HeffX2K2 = _chainJVals[LHS_spin] * mx2LHS + _chainJVals[spin] * mx2RHS + _biasFieldDriving*cos(_drivingAngFreq * t0HalfStep);
            } else {
                HeffX2K2 = _chainJVals[LHS_spin] * mx2LHS + _chainJVals[spin] * mx2RHS;
            }

            HeffY2K2 = _chainJVals[LHS_spin] * my2LHS + _chainJVals[spin] * my2RHS;
            HeffZ2K2 = _chainJVals[LHS_spin] * mz2LHS + _chainJVals[spin] * mz2RHS + GV.GetBiasField();

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

void Numerical_Methods_Class::RK2Midpoint() {

    progressbar bar(100);

    InformUserOfCodeType();

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    //std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    std::ofstream SystemAllValuesPart1(GV.GetFilePath() + "rk2_analysisFile1_" + GV.GetFileNameBase() + ".csv");
    std::ofstream SystemAllValuesPart2(GV.GetFilePath() + "rk2_analysisFile2_" + GV.GetFileNameBase() + ".csv");
    //CreateFileHeader(mxRK2File, _saveAllSpins, _onlyShowFinalState);

    //std::ofstream myRK2File(GV.GetFilePath()+"rk2_my_"+GV.GetFileNameBase()+".csv");
    //CreateFileHeader(myRK2File, _saveAllSpins, _onlyShowFinalState);
    //std::ofstream mzRK2File(GV.GetFilePath()+"rk2_mz_"+GV.GetFileNameBase()+".csv");
    //CreateFileHeader(mzRK2File, _saveAllSpins, _onlyShowFinalState);

    std::ofstream RK2MidPointMAnalysis(GV.GetFilePath() + "rk2_mx_analysis_" + GV.GetFileNameBase() + ".csv");
    RK2MidPointMAnalysis << "Iteration" << ", " << "Max M Value" << ", " << "Min M Value" << ", " << "Av M Value" << ", " << "No. over 1.0" << ", " << "Site Listing" << std::endl;

    SystemAllValuesPart1 << "Spin,Iteration,t0,hX0,hY0,hZ0,mxK1,myK1,mzK1,mx1,my1,mz1\n";
    SystemAllValuesPart2 << "Spin,Iteration,t0h,hX1, hY1, hZ1,mxK2,myK2,mzK2,mx2,my2,mz2,M Norm, Delta Norm\n";

    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {
        
        if (iteration % (_iterationEnd / 100) == 0)
            bar.update(); // Doesn't work on Windows due to different compiler

        double maxMSum = 1.0, minMSum = 1.0, averageMSum = 0.0;
        int counter = 0;
        std::vector<int> sitesLargerThanMSum{0};

        SetShockwaveConditions(iteration);

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

            // The m-components for the first stage for the: current spin site (CUR); site to the left (LHS); site to the right (RHS)
            double mx0MID = _mx0[spin], mx0LHS = _mx0[spinLHS], mx0RHS = _mx0[spinRHS];
            double my0MID = _my0[spin], my0LHS = _my0[spinLHS], my0RHS = _my0[spinRHS];
            double mz0MID = _mz0[spin], mz0LHS = _mz0[spinLHS], mz0RHS = _mz0[spinRHS];

            double hX0; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) 
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                hX0 = _chainJVals[spinLHS] * mx0LHS + _chainJVals[spin] * mx0RHS + _biasFieldDriving * cos(_drivingAngFreq * t0);
            else 
                // All spins along x which are not within the driving region
                hX0 = _chainJVals[spinLHS] * mx0LHS + _chainJVals[spin] * mx0RHS;

            // No changes are made to the effective field in the y-direction
            double hY0 = _chainJVals[spinLHS] * my0LHS + _chainJVals[spin] * my0RHS;
            // The static bias field is applied in the z-direction
            double hZ0 = _chainJVals[spinLHS] * mz0LHS + _chainJVals[spin] * mz0RHS + GV.GetBiasField();

            double mxK1, myK1, mzK1; // These are the estimations of the slopes at the beginning of the interval
            if (_useLLG) {
                /**
                 * The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters
                 * for the first stage of RK2.
                 */
                mxK1 = _gyroMagConst * (- (_gilbertVector[spin] * hY0 * mx0MID * my0MID) + hY0 * mz0MID - hZ0 * (my0MID + _gilbertVector[spin]*mx0MID*mz0MID) + _gilbertVector[spin] * hX0 * (pow(my0MID,2) + pow(mz0MID,2)));
                myK1 = _gyroMagConst * (-(hX0*mz0MID) + hZ0 * (mx0MID - _gilbertVector[spin] * my0MID * mz0MID) + _gilbertVector[spin] * (hY0 * pow(mx0MID,2) - hX0 * mx0MID * my0MID + hY0 * pow(mz0MID,2)));
                mzK1 = _gyroMagConst * (hX0 * my0MID + _gilbertVector[spin] * hZ0*(pow(mx0MID,2) + pow(my0MID,2)) - _gilbertVector[spin]*hX0*mx0MID*mz0MID - hY0 * (mx0MID + _gilbertVector[spin] * my0MID * mz0MID));
            } else {
                /**
                 * The magnetic moment components' coupled equations (obtained from the torque equation) with the
                 * parameters for the first stage of RK2.
                 * */
                mxK1 = -1 * _gyroMagConst * (my0MID * hZ0 - mz0MID * hY0);
                myK1 = _gyroMagConst * (mx0MID * hZ0 - mz0MID * hX0);
                mzK1 = -1 * _gyroMagConst * (mx0MID * hY0 - my0MID * hX0);
            }

            // Find (m0 + k1/2) for each spin, which is used in the next stage.
            mx1[spin] = mx0MID + mxK1 * _stepsizeHalf;
            my1[spin] = my0MID + myK1 * _stepsizeHalf;
            mz1[spin] = mz0MID + mzK1 * _stepsizeHalf;

            SystemAllValuesPart1 << spin << ", " << iteration << ", " << _totalTime << ", " << hX0 << ", " << hY0 << ", " << hZ0 << ", " << mxK1 << ", " << myK1 << ", " << mzK1 << ", " << mx1[spin] << ", " << my1[spin] << ", " << mz1[spin] << "\n";
        }
        SystemAllValuesPart1 << "\n\n\n";
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
            double mx1CUR = mx1[spin], mx2LHS = mx1[spinLHS], mx2RHS = mx1[spinRHS];
            double my1CUR = my1[spin], my2LHS = my1[spinLHS], my2RHS = my1[spinRHS];
            double mz1CUR = mz1[spin], mz2LHS = mz1[spinLHS], mz2RHS = mz1[spinRHS];
            
            double hX1; // The effective field (H_eff) component acting upon each spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS)
                // If a spin is driven during Stage 1 of an RK method, then it must be driven throughout the rest of the method's stages. Note the different time value used
                hX1 = _chainJVals[spinLHS] * mx2LHS + _chainJVals[spin] * mx2RHS + _biasFieldDriving * cos(_drivingAngFreq * t0HalfStep);
            else
                hX1 = _chainJVals[spinLHS] * mx2LHS + _chainJVals[spin] * mx2RHS;

            double hY1 = _chainJVals[spinLHS] * my2LHS + _chainJVals[spin] * my2RHS;
            double hZ1 = _chainJVals[spinLHS] * mz2LHS + _chainJVals[spin] * mz2RHS + GV.GetBiasField();

            double mxK2, myK2, mzK2;
            if (_useLLG) {
                // The magnetic moment components' coupled equations (obtained from LLG equation)
                mxK2 = _gyroMagConst * (- (_gilbertVector[spin] * hY1 * mx1CUR * my1CUR) + hY1 * mz1CUR - hZ1 * (my1CUR + _gilbertVector[spin]*mx1CUR*mz1CUR) + _gilbertVector[spin] * hX1 * (pow(my1CUR,2) + pow(mz1CUR,2)));
                myK2 = _gyroMagConst * (-(hX1*mz1CUR) + hZ1 * (mx1CUR - _gilbertVector[spin] * my1CUR * mz1CUR) + _gilbertVector[spin] * (hY1 * pow(mx1CUR,2) - hX1 * mx1CUR * my1CUR + hY1 * pow(mz1CUR,2)));
                mzK2 = _gyroMagConst * (hX1 * my1CUR + _gilbertVector[spin] * hZ1*(pow(mx1CUR,2) + pow(my1CUR,2)) - _gilbertVector[spin]*hX1*mx1CUR*mz1CUR - hY1 * (mx1CUR + _gilbertVector[spin] * my1CUR * mz1CUR));
            } else {
                // The magnetic moment components' coupled equations (obtained from the torque equation)
                mxK2 = -1 * _gyroMagConst * (my1CUR * hZ1 - mz1CUR * hY1);
                myK2 =      _gyroMagConst * (mx1CUR * hZ1 - mz1CUR * hX1);
                mzK2 = -1 * _gyroMagConst * (mx1CUR * hY1 - my1CUR * hX1);
            }

            mx2[spin] = _mx0[spin] + mxK2 * _stepsize;
            my2[spin] = _my0[spin] + myK2 * _stepsize;
            mz2[spin] = _mz0[spin] + mzK2 * _stepsize;

            double mSumTotal = mx2[spin] + my2[spin] + mz2[spin];
            SystemAllValuesPart2 << spin << ", " << iteration << ", " << t0HalfStep << ", " << hX1 << ", " << hY1 << ", " << hZ1 << ", " << mxK2 << ", " << myK2 << ", " << mzK2 << ", " << mx2[spin] << ", " << my2[spin] << ", " << mz2[spin] << ", " << mSumTotal << ", " << 1.0 - mSumTotal << "\n";
            // std::cout << mz2[spin] << " | " << my2[spin] << " | " << mz2[spin] << std::endl;

            if (_shouldTrackMValues) {
                if (spin > _numSpinsDamped && spin <= _numSpinsDamped + _numSpinsInChain) {
                    double mSumTotal = mx2[spin] + my2[spin] + mz2[spin];

                    if (mSumTotal > maxMSum) {
                        maxMSum = mSumTotal;
                        _maxMSummation = mSumTotal;
                        sitesLargerThanMSum.push_back(spin);
                        counter++;
                    }

                    if (mSumTotal < minMSum) {
                        minMSum = mSumTotal;
                    }
                    averageMSum += mSumTotal;
                }

                double mSumNorm = sqrt(pow(mx2[spin], 2) + pow(my2[spin], 2) + pow(mz2[spin], 2));

                mx2[spin] /= mSumNorm;
                my2[spin] /= mSumNorm;
                mz2[spin] /= mSumNorm;
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.
        SystemAllValuesPart2 << "\n\n\n";
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

        //SaveDataToFile(_saveAllSpins, mxRK2File, mx2, iteration, _onlyShowFinalState);
        //SaveDataToFile(_saveAllSpins, myRK2File, my2, iteration, _onlyShowFinalState);
        //SaveDataToFile(_saveAllSpins, mzRK2File, mz2, iteration, _onlyShowFinalState);

        if (_shouldTrackMValues) {
            if (iteration % static_cast<int>(_iterationEnd / _numberOfDataPoints) == 0) {
                RK2MidPointMAnalysis << iteration << ", " << maxMSum << ", " << minMSum << ", " << (averageMSum) / _numSpinsInChain << ", " << counter << ", ";
                for (int i: sitesLargerThanMSum) { RK2MidPointMAnalysis << i << ", "; }
                RK2MidPointMAnalysis << std::endl;
            }
        }

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        _mx0 = mx2;
        _my0 = my2;
        _mz0 = mz2;

        if (iteration == 10)
            exit(0);
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    //mxRK2File.close();
    //myRK2File.close();
    //mzRK2File.close();
    RK2MidPointMAnalysis.close();
    SystemAllValuesPart1.close();
    SystemAllValuesPart2.close();

    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "\nmax value is: " << _maxMSummation << std::endl;
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::RK2LLGTestbed() {

    progressbar bar(100);

    double maxVal = 1.0;
    // Notifies the user of what code they are running.
    InformUserOfCodeType();

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath()+"rk2_mx_"+GV.GetFileNameBase()+".csv");
    //std::ofstream myRK2File(GV.GetFilePath()+"rk2_my_"+GV.GetFileNameBase()+".csv");
    //std::ofstream mzRK2File(GV.GetFilePath()+"rk2_mz_"+GV.GetFileNameBase()+".csv");

    std::ofstream mxRK2Analysis(GV.GetFilePath()+"rk2_mx_analysis_"+GV.GetFileNameBase()+".csv");
    mxRK2Analysis << "Iteration" << ", " << "Max M Value" << ", " << "Min M Value" << ", " << "Av M Value" << ", " << "No. over 1.0" << ", " << "Site Listing" << std::endl;

    CreateFileHeader(mxRK2File, _saveAllSpins, _onlyShowFinalState);
    //CreateFileHeader(myRK2File, _saveAllSpins, _onlyShowFinalState);
    //CreateFileHeader(mzRK2File, _saveAllSpins, _onlyShowFinalState);

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or half-step (h/2)*/
    for (int iteration = _iterationStart; iteration <= _iterationEnd; iteration++) {
        // Doesn't work on Windows
        if (iteration % (_iterationEnd / 100) == 0) {
            bar.update();
        }

        double maxMVal = 1.0, minMVal = 1.0, avMVal = 0.0;
        int counter = 0;
        std::vector<int> largesites {0};

        SetShockwaveConditions(iteration);

        _totalTime += _stepsize;
        double t0 = _totalTime; // The initial time of the iteration, and the time at the first stage; the start of the interval and the first step of RK2. Often called 't0' in literature
        double t0HalfStep = _totalTime + _stepsizeHalf; // The time at the midpoint of the interval; where the second stage of RK2 is. Often called 't0*h/2' in literature

        // The estimate of the slope for the x-axis magnetic moment component at the midpoint
        std::vector<double> mxEstMid(GV.GetNumSpins() + 2, 0);
        // The estimate of the slope for the y-axis magnetic moment component at the midpoint
        std::vector<double> myEstMid(GV.GetNumSpins() + 2, 0);
        // The estimate of the slope for the z-axis magnetic moment component at the midpoint
        std::vector<double> mzEstMid(GV.GetNumSpins() + 2, 0);

        // Loop the 0th and final spins as they will always be zero-valued
        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            /* The first stage is based upon finding the value of the slope at the beginning of the interval (k1). This
             * stage takes the start conditions as an input, and substitutes them into the LLG equation. */
            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The mx components for the first stage for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mx1 = _mx0[spin], mx1LHS = _mx0[spinLHS], mx1RHS = _mx0[spinRHS];
            // The my components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double my1 = _my0[spin], my1LHS = _my0[spinLHS], my1RHS = _my0[spinRHS];
            // The mz components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mz1 = _mz0[spin], mz1LHS = _mz0[spinLHS], mz1RHS = _mz0[spinRHS];

            double mxK1, myK1, mzK1; // These are the estimations of the slopes at the beginning of the interval for each magnetic moment component
            double heffXK1, heffYK1, heffZK1; // The effective field component acting upon each spin

            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                heffXK1 = _chainJVals[spinLHS] * mx1LHS + _chainJVals[spin] * mx1RHS + _biasFieldDriving * cos(_drivingAngFreq * t0);
            } else {
                // The ELSE IF statement includes all spins along x which are not within the driving region
                heffXK1 = _chainJVals[spinLHS] * mx1LHS + _chainJVals[spin] * mx1RHS;
            }

            // No changes are made to the effective field in the y-direction
            heffYK1 = _chainJVals[spinLHS] * my1LHS + _chainJVals[spin] * my1RHS;
            // The bias field is applied in the z-direction and so it contributes to the effective field in the z-direction
            heffZK1 = _chainJVals[spinLHS] * mz1LHS + _chainJVals[spin] * mz1RHS + GV.GetBiasField();

            if (_useLLG) {
                /* The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters
                 * for the first stage of RK2. */
                mxK1 = _gyroMagConst * (- (_gilbertVector[spin] * heffYK1 * mx1 * my1) + heffYK1 * mz1 - heffZK1 * (my1 + _gilbertVector[spin]*mx1*mz1) + _gilbertVector[spin] * heffXK1 * (pow(my1,2) + pow(mz1,2)));
                myK1 = _gyroMagConst * (-(heffXK1*mz1) + heffZK1 * (mx1 - _gilbertVector[spin] * my1 * mz1) + _gilbertVector[spin] * (heffYK1 * pow(mx1,2) - heffXK1 * mx1 * my1 + heffYK1 * pow(mz1,2)));
                mzK1 = _gyroMagConst * (heffXK1 * my1 + _gilbertVector[spin] * heffZK1*(pow(mx1,2) + pow(my1,2)) - _gilbertVector[spin]*heffXK1*mx1*mz1 - heffYK1 * (mx1 + _gilbertVector[spin] * my1 * mz1));
            } else {
                /* The magnetic moment components' coupled equations (obtained from the torque equation) with the
                 * parameters for the first stage of RK2. */
                mxK1 = -1 * _gyroMagConst * (my1 * heffZK1 - mz1 * heffYK1);
                myK1 = _gyroMagConst * (mx1 * heffZK1 - mz1 * heffXK1);
                mzK1 = -1 * _gyroMagConst * (mx1 * heffYK1 - my1 * heffXK1);
            }

            mxEstMid[spin] = mx1 + mxK1 * _stepsizeHalf;
            myEstMid[spin] = my1 + myK1 * _stepsizeHalf;
            mzEstMid[spin] = mz1 + mzK1 * _stepsizeHalf;
        }

        // The estimate of the mx value for the next iteration of iteration calculated using the RK2 Midpoint rule
        std::vector<double> mxNextVal(GV.GetNumSpins() + 2,0);
        // The estimate of the mY value for the next iteration of iteration
        std::vector<double> myNextVal(GV.GetNumSpins() + 2,0);
        // The estimate of the mz value for the next iteration of iteration
        std::vector<double> mzNextVal(GV.GetNumSpins() + 2,0);

        for (int spin = 1; spin <= GV.GetNumSpins(); spin++) {
            /* The second stage uses the previously found k1 value, as well as the initial conditions, to determine the
             * value of the slope (k2) at the midpoint. In RK2, the values of k1 and k2 can then be jointly used to
             * estimate the next point of the function through a weighted average of k1 & k2.
             *
             * In this loop the definitions of variables follow a similar format to Stage1.*/

            int spinLHS = spin - 1, spinRHS = spin + 1;
            double mx2 = mxEstMid[spin], mx2LHS = mxEstMid[spinLHS], mx2RHS = mxEstMid[spinRHS];
            double my2 = myEstMid[spin], my2LHS = myEstMid[spinLHS], my2RHS = myEstMid[spinRHS];
            double mz2 = mzEstMid[spin], mz2LHS = mzEstMid[spinLHS], mz2RHS = mzEstMid[spinRHS];

            double mxK2, myK2, mzK2;
            double heffXK2, heffYK2, heffZK2;

            // Driving region must be consistently applied at every stage of the RK2 method
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS)
                // If a spin is driven during Stage 1 of an RK method, then it must be driven throughout the rest of the method's stages
                heffXK2 = _chainJVals[spinLHS] * mx2LHS + _chainJVals[spin] * mx2RHS + _biasFieldDriving * cos(_drivingAngFreq * t0HalfStep);
            else
                heffXK2 = _chainJVals[spinLHS] * mx2LHS + _chainJVals[spin] * mx2RHS;

            heffYK2 = _chainJVals[spinLHS] * my2LHS + _chainJVals[spin] * my2RHS;
            heffZK2 = _chainJVals[spinLHS] * mz2LHS + _chainJVals[spin] * mz2RHS + GV.GetBiasField();

            if (_useLLG) {
                /* The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters
                 * for the second stage of RK2. */
                mxK2 = _gyroMagConst * (- (_gilbertVector[spin] * heffYK2 * mx2 * my2) + heffYK2 * mz2 - heffZK2 * (my2 + _gilbertVector[spin]*mx2*mz2) + _gilbertVector[spin] * heffXK2 * (pow(my2,2) + pow(mz2,2)));
                myK2 = _gyroMagConst * (-(heffXK2*mz2) + heffZK2 * (mx2 - _gilbertVector[spin] * my2 * mz2) + _gilbertVector[spin] * (heffYK2 * pow(mx2,2) - heffXK2 * mx2 * my2 + heffYK2 * pow(mz2,2)));
                mzK2 = _gyroMagConst * (heffXK2 * my2 + _gilbertVector[spin] * heffZK2*(pow(mx2,2) + pow(my2,2)) - _gilbertVector[spin]*heffXK2*mx2*mz2 - heffYK2 * (mx2 + _gilbertVector[spin] * my2 * mz2));
            } else {
                /* The magnetic moment components' coupled equations (obtained from the torque equation) with the
                 * parameters for the second stage of RK2. */
                mxK2 = -1 * _gyroMagConst * (my2 * heffZK2 - mz2 * heffYK2);
                myK2 =      _gyroMagConst * (mx2 * heffZK2 - mz2 * heffXK2);
                mzK2 = -1 * _gyroMagConst * (mx2 * heffYK2 - my2 * heffXK2);
            }

            mxNextVal[spin] = _mx0[spin] + mxK2 * _stepsize;
            myNextVal[spin] = _my0[spin] + myK2 * _stepsize;
            mzNextVal[spin] = _mz0[spin] + mzK2 * _stepsize;

            if (spin > _numSpinsDamped && spin <= _numSpinsDamped + _numSpinsInChain) {
                double mValsTotal = mxNextVal[spin] + myNextVal[spin] + mzNextVal[spin];
                // if (mValsTotal > 1.1) {std::cout << "Iter: " << iteration << " | Spin: " << spin << " | MVal: " << mValsTotal << std::endl;}
                if (mValsTotal > maxMVal) {
                    maxMVal = mValsTotal;
                    maxVal = mValsTotal;
                    largesites.push_back(spin);
                    counter++;
                }

                if (mValsTotal < minMVal) {
                    minMVal = mValsTotal;
                }
                avMVal += mValsTotal;
            }

            double mValsNorm = sqrt( pow(mxNextVal[spin], 2) + pow(myNextVal[spin], 2) + pow(mzNextVal[spin], 2));

            mxNextVal[spin] /= mValsNorm;
            myNextVal[spin] /= mValsNorm;
            mzNextVal[spin] /= mValsNorm;

        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.


        // Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
        // these between loop iterations sometimes led to incorrect values cropping up.
        _mx0.clear();
        _my0.clear();
        _mz0.clear();
        mxEstMid.clear();
        myEstMid.clear();
        mzEstMid.clear();

        SaveDataToFile(_saveAllSpins, mxRK2File, mxNextVal, iteration, _onlyShowFinalState);
        //SaveDataToFile(_saveAllSpins, myRK2File, myNextVal, iteration, _onlyShowFinalState);
        //SaveDataToFile(_saveAllSpins, mzRK2File, mzNextVal, iteration, _onlyShowFinalState);
        if (iteration % static_cast<int>(_iterationEnd / _numberOfDataPoints) == 0) {
            mxRK2Analysis << iteration << ", " << maxMVal << ", " << minMVal << ", " << (avMVal) / _numSpinsInChain << ", " << counter << ", ";
            for (int i: largesites) {
                if (i > _numSpinsDamped && i <= _numSpinsDamped + _numSpinsInChain){mxRK2Analysis << i << ", ";}
            }
            mxRK2Analysis << std::endl;
        }
        /* Sets the final value of the current iteration of the loop (y_(n+1) in textbook's notation) to be the starting
         * value of the next iteration (y_n) */
        _mx0 = mxNextVal;
        _my0 = myNextVal;
        _mz0 = mzNextVal;
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();
    mxRK2Analysis.close();
    //myRK2File.close();
    //mzRK2File.close();
    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "\nmax value is: " << maxVal << std::endl;
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::RK4() {
    double appliedStaticField = 0.1;  // [T]
    double appliedDynamicField = 3E-3;  // [T]
    int numberOfSpins = 4000;

    double exchangeMin = 43.5, exchangeMax = 132.0; // [T]
    double drivingFrequency = 42.5 * 1e9, angularDrivingFrequency = 2.0 * M_PI * drivingFrequency;
    double gamma = 29.2E9 * (2 * M_PI); // [Hz / T]

    int drivingRegionLHS = 1, drivingRegionWidth = 200, drivingRegionRHS = drivingRegionLHS + drivingRegionWidth;

    double totalTime = 0;
    double stepsize = 1e-15, halfStepsize = stepsize / 2.0;
    double maxIterations = 7e5;
    int numberOfDatapoints = 100;
    int numberOfPairs = numberOfSpins - 1;

    double maxMVal = 1.0;

    //########################################################################################################################
    double mxInit = 0.0, myInit = 0.0, mzInit = 1.0;
    std::vector<double> mxInitCond(numberOfSpins, mxInit), myInitCond(numberOfSpins, myInit), mzInitCond(numberOfSpins, mzInit);

    std::vector<double> mxStart{0}, myStart{0}, mzStart{0};

    // Appends initial conditions to the vectors
    mxStart.insert(mxStart.end(), mxInitCond.begin(), mxInitCond.end());
    myStart.insert(myStart.end(), myInitCond.begin(), myInitCond.end());
    mzStart.insert(mzStart.end(), mzInitCond.begin(), mzInitCond.end());

    // This zero is the (N+1)th spin on the RHS of the chain
    mxStart.push_back(0);
    myStart.push_back(0);
    mzStart.push_back(0);

    // The linearly spaced vector is saved as the class member '_chainJVals' simply to increase code readability
    std::vector<double> exchangeArray = _chainJVals;
    //########################################################################################################################

    std::string filepath = "/Users/cameronmceleney/CLionProjects/Data/2022-05-18/";
    std::ofstream mxRK2File(filepath + "rk4_mx_" + std::to_string(numberOfSpins) + ".csv");

    for (int iteration = 0; iteration <= maxIterations; iteration++) {
        totalTime += stepsize;
        double t0 = totalTime, t0Half = totalTime + halfStepsize, t0h = totalTime + stepsize;

        std::vector<double> mxEstStep1(numberOfSpins+2, 0), myEstStep1(numberOfSpins+2, 0), mzEstStep1(numberOfSpins+2, 0);
        std::vector<double> mxk1Est (numberOfSpins + 2, 0), myk1Est (numberOfSpins + 2, 0), mzk1Est (numberOfSpins + 2, 0);

        for (int spin = 1; spin <= numberOfSpins; spin++) {
            // RK4 Step 1
            int spinLHS = spin - 1, spinRHS = spin + 1;

            double mx1 = mxStart[spin], mx1LHS = mxStart[spinLHS], mx1RHS = mxStart[spinRHS];
            double my1 = myStart[spin], my1LHS = myStart[spinLHS], my1RHS = myStart[spinRHS];
            double mz1 = mzStart[spin], mz1LHS = mzStart[spinLHS], mz1RHS = mzStart[spinRHS];

            double hX1, hY1, hZ1;

            if (spin >= drivingRegionLHS && spin <= drivingRegionRHS) {
                hX1 = exchangeArray[spinLHS] * mx1LHS + exchangeArray[spin] * mx1RHS + appliedDynamicField * cos(angularDrivingFrequency * t0);
            } else {
                hX1 = exchangeArray[spinLHS] * mx1LHS + exchangeArray[spin] * mx1RHS;
            }
            hY1 = exchangeArray[spinLHS] * my1LHS + exchangeArray[spin] * my1RHS;
            hZ1 = exchangeArray[spinLHS] * mz1LHS + exchangeArray[spin] * mz1RHS + appliedStaticField;

            double mxK1 = -1.0 * gamma * (my1 * hZ1 - mz1 * hY1);
            double myK1 = +1.0 * gamma * (mx1 * hZ1 - mz1 * hX1);
            double mzK1 = -1.0 * gamma * (mx1 * hY1 - my1 * hX1);

            mxk1Est[spin] = mxK1;
            myk1Est[spin] = myK1;
            mzk1Est[spin] = mzK1;

            mxEstStep1[spin] = mx1 + mxK1 * halfStepsize;
            myEstStep1[spin] = my1 + myK1 * halfStepsize;
            mzEstStep1[spin] = mz1 + mzK1 * halfStepsize;
        }

        std::vector<double> mxEstStep2 (numberOfSpins + 2, 0), myEstStep2 (numberOfSpins + 2, 0), mzEstStep2 (numberOfSpins + 2, 0);
        std::vector<double> mxk2Est (numberOfSpins + 2, 0), myk2Est (numberOfSpins + 2, 0), mzk2Est (numberOfSpins + 2, 0);

        for (int spin = 1; spin <= numberOfSpins; spin++) {
            // RK4 Step 2
            int spinLHS = spin - 1, spinRHS = spin + 1;

            double mx2 = mxEstStep1[spin], mx2LHS = mxEstStep1[spinLHS], mx2RHS = mxEstStep1[spinRHS];
            double my2 = myEstStep1[spin], my2LHS = myEstStep1[spinLHS], my2RHS = myEstStep1[spinRHS];
            double mz2 = mzEstStep1[spin], mz2LHS = mzEstStep1[spinLHS], mz2RHS = mzEstStep1[spinRHS];

            double hX2, hY2, hZ2;
            if (spin >= drivingRegionLHS && spin <= drivingRegionRHS) {
                hX2 = exchangeArray[spinLHS] * mx2LHS + exchangeArray[spin] * mx2RHS + appliedDynamicField * cos(angularDrivingFrequency * t0Half);
            } else {
                hX2 = exchangeArray[spinLHS] * mx2LHS + exchangeArray[spin] * mx2RHS;
            }
            hY2 = exchangeArray[spinLHS] * my2LHS + exchangeArray[spin] * my2RHS;
            hZ2 = exchangeArray[spinLHS] * mz2LHS + exchangeArray[spin] * mz2RHS + appliedStaticField;

            double mxK2 = -1.0 * gamma * (my2 * hZ2 - mz2 * hY2);
            double myK2 = +1.0 * gamma * (mx2 * hZ2 - mz2 * hX2);
            double mzK2 = -1.0 * gamma * (mx2 * hY2 - my2 * hX2);

            mxk2Est[spin] = mxK2;
            myk2Est[spin] = myK2;
            mzk2Est[spin] = mzK2;

            mxEstStep2[spin] = mxStart[spin] + mxK2 * halfStepsize;
            myEstStep2[spin] = myStart[spin] + myK2 * halfStepsize;
            mzEstStep2[spin] = mzStart[spin] + mzK2 * halfStepsize;
        }

        std::vector<double> mxEstStep3 (numberOfSpins + 2, 0), myEstStep3 (numberOfSpins + 2, 0), mzEstStep3 (numberOfSpins + 2, 0);
        std::vector<double> mxk3Est (numberOfSpins + 2, 0), myk3Est (numberOfSpins + 2, 0), mzk3Est (numberOfSpins + 2, 0);

        for (int spin = 1; spin <= numberOfSpins; spin++) {
            // RK4 Step 3
            int spinLHS = spin - 1, spinRHS = spin + 1;

            double mx3 = mxEstStep2[spin], mx3LHS = mxEstStep2[spinLHS], mx3RHS = mxEstStep2[spinRHS];
            double my3 = myEstStep2[spin], my3LHS = myEstStep2[spinLHS], my3RHS = myEstStep2[spinRHS];
            double mz3 = mzEstStep2[spin], mz3LHS = mzEstStep2[spinLHS], mz3RHS = mzEstStep2[spinRHS];

            double hX3, hY3, hZ3;
            if (spin >= drivingRegionLHS && spin <= drivingRegionRHS) {
                hX3 = exchangeArray[spinLHS] * mx3LHS + exchangeArray[spin] * mx3RHS + appliedDynamicField * cos(angularDrivingFrequency * t0Half);
            } else {
                hX3 = exchangeArray[spinLHS] * mx3LHS + exchangeArray[spin] * mx3RHS;
            }
            hY3 = exchangeArray[spinLHS] * my3LHS + exchangeArray[spin] * my3RHS;
            hZ3 = exchangeArray[spinLHS] * mz3LHS + exchangeArray[spin] * mz3RHS + appliedStaticField;

            double mxK3 = -1.0 * gamma * (my3 * hZ3 - mz3 * hY3);
            double myK3 = +1.0 * gamma * (mx3 * hZ3 - mz3 * hX3);
            double mzK3 = -1.0 * gamma * (mx3 * hY3 - my3 * hX3);

            mxk3Est[spin] = mxK3;
            myk3Est[spin] = myK3;
            mzk3Est[spin] = mzK3;

            mxEstStep3[spin] = mxStart[spin] + mxK3 * stepsize;
            myEstStep3[spin] = myStart[spin] + myK3 * stepsize;
            mzEstStep3[spin] = mzStart[spin] + mzK3 * stepsize;
        }

        std::vector<double> mxEstFinal (numberOfSpins + 2, 0), myEstFinal (numberOfSpins + 2, 0), mzEstFinal (numberOfSpins + 2, 0);
        for (int spin = 1; spin <= numberOfSpins; spin++) {
            // RK4 Step 4
            int spinLHS = spin - 1, spinRHS = spin + 1;

            double mx4 = mxEstStep3[spin], mx4LHS = mxEstStep3[spinLHS], mx4RHS = mxEstStep3[spinRHS];
            double my4 = myEstStep3[spin], my4LHS = myEstStep3[spinLHS], my4RHS = myEstStep3[spinRHS];
            double mz4 = mzEstStep3[spin], mz4LHS = mzEstStep3[spinLHS], mz4RHS = mzEstStep3[spinRHS];

            double hX4, hY4, hZ4;

            if (spin >= drivingRegionLHS && spin <= drivingRegionRHS) {
                hX4 = exchangeArray[spinLHS] * mx4LHS + exchangeArray[spin] * mx4RHS + appliedDynamicField * cos(angularDrivingFrequency * t0h);
            } else {
                hX4 = exchangeArray[spinLHS] * mx4LHS + exchangeArray[spin] * mx4RHS;
            }
            hY4 = exchangeArray[spinLHS] * my4LHS + exchangeArray[spin] * my4RHS;
            hZ4 = exchangeArray[spinLHS] * mz4LHS + exchangeArray[spin] * mz4RHS + appliedStaticField;

            double mxK4 = -1.0 * gamma * (my4 * hZ4 - mz4 * hY4);
            double myK4 = +1.0 * gamma * (mx4 * hZ4 - mz4 * hX4);
            double mzK4 = -1.0 * gamma * (mx4 * hY4 - my4 * hX4);

            mxEstFinal[spin] = mxStart[spin] + (mxk2Est[spin]  + 2.0*mxk2Est[spin] + 2.0*mxk2Est[spin] + mxK4) * (stepsize / 6.0);
            myEstFinal[spin] = myStart[spin] + (myk2Est[spin] + 2.0*myk2Est[spin] + 2.0*myk2Est[spin] + myK4) * (stepsize / 6.0);
            mzEstFinal[spin] = mzStart[spin] + (mzk2Est[spin] + 2.0*mzk2Est[spin] + 2.0*mzk2Est[spin] + mzK4) * (stepsize / 6.0);


            double mValsTotal = mxEstFinal[spin] + myEstFinal[spin] + mzEstFinal[spin];
            if (mValsTotal > maxMVal) {
                maxMVal = mValsTotal;
            }
        }

        if (iteration % static_cast<int>(maxIterations / numberOfDatapoints) == 0) {
            //if (iteration == _iterationEnd) {
            for (int i = 0; i <= numberOfSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    // Print current time
                    mxRK2File << (iteration * stepsize) << ",";

                else if (i == numberOfSpins)
                    // Ensures that the final line doesn't contain a comma.
                    mxRK2File << mxEstStep2[i] << std::flush;

                else
                    // For non-special values, write the data.
                    mxRK2File << mxEstStep2[i] << ",";
            }
            // Take new line after current row is finished being written.
            mxRK2File << std::endl;
        }


        mxStart = mxEstFinal;
        myStart = myEstFinal;
        mzStart = mzEstFinal;
    }

    mxRK2File.close();

    std::cout << "Max M Value: " << maxMVal << std::endl;
}

void Numerical_Methods_Class::InformUserOfCodeType() {
    /**
     * Informs the user of the code type they are running, including: solver type; special modules.
     */
    if (_useLLG)
        std::cout << "\nYou are running the RK2 Spinchains (LLG) code";
    else
        std::cout << "\nYou are running the RK2 Spinchains (Torque) code";

    if (_hasShockwave)
        std::cout << " with shockwave module.\n";
    else
        std::cout << ".\n";

}
void Numerical_Methods_Class::CreateFileHeader(std::ofstream &outputFileName, bool &areAllSpinBeingSaved, bool &onlyShowFinalState) {
    /**
     * Write all non-data information to the output file.
     */
    outputFileName << "Key Data\n";

    outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

    outputFileName << "Using LLG," << _useLLG << ", " << "Using Shockwave," << _hasShockwave << ", " << "Drive from LHS," << _lhsDrive << "\n";

    outputFileName << "\n";

    outputFileName << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                      "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site, Driving Region Width,"
                      "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                      "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, Stepsize (h),Gilbert Damping Factor, Gyromagnetic Ratio (2Pi*Y),"
                      "Shockwave Gradient Time (s), Shockwave Application Time (s)"
                      "\n";

    outputFileName << GV.GetBiasField() << ", " << _biasFieldDriving << ", " << _shockwaveScaling << ", " << _biasFieldDriving * _shockwaveScaling << ", "
                   << _drivingFreq << ", " << _drivingRegionLHS - _numSpinsDamped << ", " << _drivingRegionRHS - _numSpinsDamped << ", " << _drivingRegionWidth << ", "
                   << _maxSimTime << ", " << GV.GetExchangeMinVal() << ", " << GV.GetExchangeMaxVal() << ", " << _iterationEnd << ", " <<  _numberOfDataPoints << ", "
                   << _numSpinsInChain << ", "  << _numSpinsDamped << ", " << _numSpinsInChain + 2 * _numSpinsDamped << ", " << _stepsize << ", " << _gilbertConst << ", " << _gyroMagConst << ", "
                   << _iterToBeginShockwave << ", " << _shockwaveGradientTime * _stepsize
                   << "\n";

    outputFileName << "\n";

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments );
    outputFileName << "Note(s):," << notesComments << "\n"; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n";

    outputFileName << "\n";

    CreateColumnHeaders(outputFileName, areAllSpinBeingSaved, onlyShowFinalState);

    std::cout << "\n";
}
void Numerical_Methods_Class::CreateColumnHeaders(std::ofstream &outputFileName, bool &areAllSpinBeingSaved, bool &onlyShowFinalState) {
    /**
     * Creates the column headers for each spin site simulated. This code can change often, so compartmentalising it in
     * a separate function is necessary to reduce bugs.
     */
    if (areAllSpinBeingSaved or onlyShowFinalState) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s]";
        for (int i = 1; i <= GV.GetNumSpins(); i++) {
            outputFileName << "," << i;
        }
        outputFileName << std::endl;

    } else if (_fixedPoints) {
        /*outputFileName << "Time" << ", "
                       << _drivingRegionLHS << ","
                       << static_cast<int>(_drivingRegionWidth / 2.0) << ","
                       << _drivingRegionRHS << ","
                       << static_cast<int>(1500) << ","
                       << static_cast<int>(2500) << ","
                       << static_cast<int>(3500) << ","
                       << GV.GetNumSpins() << std::endl;*/

        outputFileName << "Time" << ", "
                       << static_cast<int>(400) << ","
                       << static_cast<int>(1500) << ","
                       << static_cast<int>(3000) << ","
                       << static_cast<int>(4500) << ","
                       << static_cast<int>(5600) << std::endl;
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
void Numerical_Methods_Class::SetShockwaveConditions(double current_iteration) {

    // If function is triggered, then the applied biasFieldDriving is increased by the scale factor _shockwaveScaling
    if (_hasShockwave and not _isShockwaveOn)
    {
        if (current_iteration >= _iterationEnd * _iterToBeginShockwave)
        {
            // Shockwave begins once simulation is a certain % complete
            _isShockwaveOn = true;
        }

        return;
    }

    if (_isShockwaveOn and not _isShockwaveAtMax)
    {
        _biasFieldDriving += _shockwaveStepsize;

        if (_biasFieldDriving >= _shockwaveMax)
        {
            _biasFieldDriving = _shockwaveMax;
            _isShockwaveAtMax = true;

        }
        return;

    }

}
void Numerical_Methods_Class::SaveDataToFile(bool &areAllSpinBeingSaved, std::ofstream &outputFileName,
                                             std::vector<double> &arrayToWrite, int &iteration, bool &onlyShowFinalState) {
    if (onlyShowFinalState) {
        // iteration >= static_cast<int>(_iterationEnd / 2.0) &&
        if (iteration % (_iterationEnd / _numberOfDataPoints) == 0) {
        //if (iteration == _iterationEnd) {
            for (int i = 0; i <= GV.GetNumSpins(); i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * _stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ",";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;
        }
    } else {
        if (areAllSpinBeingSaved)
        {
            for (int i = 0; i <= GV.GetNumSpins(); i++)
            {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * _stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        }
        else
        {
            if (iteration % (_iterationEnd / _numberOfDataPoints) == 0)
            {
                if (_fixedPoints)
                {
                    /*
                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[_drivingRegionLHS] << ","
                                   << arrayToWrite[static_cast<int>(_drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[_drivingRegionRHS] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[GV.GetNumSpins()] << std::endl;
                   */
                    outputFileName << (iteration * _stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;


                }
                else
                {
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
    }
}

void Numerical_Methods_Class::DebugOptions(std::vector<double> mxNextVal, std::vector<double> myNextVal, std::vector<double> mzNextVal, int spin, long iteration) {
        if (mxNextVal[spin] >= 5000) {
            std::cout << "Error. Value of mx was greater than 5000 at spin(" << spin << "), iter("
                      << iteration << ")." << std::flush;
            std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                      << _drivingRegionLHS << "; itermax = " << _iterationEnd << "; stepSize: "
                      << _stepsize << std::endl;
            exit(3);
        }

        if (myNextVal[spin] >= 5000) {
            std::cout << "Error. Value of my was greater than 5000 at spin(" << spin << "), iter("
                      << iteration << ")." << std::flush;
            std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                      << _drivingRegionLHS << "; itermax = " << _iterationEnd << "; stepSize: "
                      << _stepsize << std::endl;
            exit(3);
        }

        if (mzNextVal[spin] >= 5000) {
            std::cout << "Error. Value of mz was greater than 5000 at spin(" << spin << "), iter("
                      << iteration << ")." << std::flush;
            std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                      << _drivingRegionLHS << "; itermax = " << _iterationEnd << "; stepSize: "
                      << _stepsize << std::endl;
            exit(3);
        }
}
