#include "Numerical_Methods_Class.h"

void Numerical_Methods_Class::NMSetup() {

    _biasFieldDriving = 1e-7;
    _drivingFreq = 42.5 * 1e9;
    _stepsize = 1e-15; // This should be at least (1 / _drivingFreq)
    _stopIterVal = static_cast<int>(1.4e6); // 2.6e5
    _undampedNumSpins = GV.GetNumSpins();

    _hasShockwave = true;
    _iterToBeginShockwave = 0.5; // Value should be between [0.0, 1.0] inclusive.
    _shockwaveScaling = 1e6;
    _shockwaveInit = _biasFieldDriving;
    _shockwaveMax = _shockwaveInit * _shockwaveScaling;
    _shockwaveIncreaseTime = 40e3; // 1 = instantaneous application. 35e3 is 35[fs] when stepsize=1e-15
    _shockwaveStepsize = (_shockwaveMax - _shockwaveInit) / _shockwaveIncreaseTime;

    _useLLG = true;
    _lhsDrive = false;

    _onlyShowFinalState = true;
    _saveAllSpins = false;
    _fixedPoints = false;

    _gilbertLower = 1e-5;
    _gilbertUpper = 1.0;
    _numGilbert = 400;
    GV.SetNumSpins(_undampedNumSpins + 2 * _numGilbert);

    _numberOfDataPoints = 100; // Set equal to _stopIterVal to save all data

    _drivingAngFreq = 2 * M_PI * _drivingFreq;
    _numberOfSpinPairs = GV.GetNumSpins() - 1;
    _stepsizeHalf = _stepsize / 2.0;
    _maxSimTime = _stepsize * _stopIterVal;

    SetDrivingRegion(_lhsDrive);

    SetupVectors(); // Calls private class method to generate vectors needed for RK methods.

    // std::cout << "gil: " << _gilbertVector.size() << ".mz: " << _mzStartVal.size() << std::endl; std::cout << "\n\n";
}
void Numerical_Methods_Class::SetDrivingRegion(bool &useLHSDrive) {

    if (useLHSDrive)
    { //Drives from the LHS, starting at _drivingRegionLHS
        _drivingRegionLHS = _numGilbert + 1; // If RHS start, then this value should be (startStart - 1) for correct offset.
        _drivingRegionWidth = 200;// static_cast<int>(_undampedNumSpins * _regionScaling);
        _drivingRegionRHS = _drivingRegionLHS + _drivingRegionWidth;
    }
    else
    { // Drives from the RHS, starting at _drivingRegionRHS
        _drivingRegionWidth = 200;// static_cast<int>(_undampedNumSpins * _regionScaling);
        _drivingRegionRHS = GV.GetNumSpins() - _numGilbert - 100 - 1;
        //_drivingRegionRHS = (_undampedNumSpins/2) +_numGilbert + (_drivingRegionWidth / 2); // use for central drive
        _drivingRegionLHS = _drivingRegionRHS - _drivingRegionWidth - 1; // The -1 is to correct the offset
    }
}
void Numerical_Methods_Class::SetupVectors() {

    SetupVectorsExchange();
    // SetupVectorsGilbert();
    GilbertVectorsBothSides();
}
void Numerical_Methods_Class::SetupVectorsExchange() {

    LinspaceClass SpinChainExchange;

    // The linearly spaced vector is saved as the class member '_chainJVals' simply to increase code readability
    SpinChainExchange.set_values(GV.GetExchangeMinVal(), GV.GetExchangeMaxVal(), _numberOfSpinPairs, true);
    SpinChainExchange.generate_array();
    _chainJVals = SpinChainExchange.build_spinchain();

    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mxInitCond(GV.GetNumSpins(), _mxInit), myInitCond(GV.GetNumSpins(), _myInit), mzInitCond(GV.GetNumSpins(), _mzInit);
    // mxInitCond[0] = _mxInit; // Only perturb initial spin

    // Appends initial conditions to the vectors
    _mxStartVal.insert(_mxStartVal.end(), mxInitCond.begin(), mxInitCond.end());
    _myStartVal.insert(_myStartVal.end(), myInitCond.begin(), myInitCond.end());
    _mzStartVal.insert(_mzStartVal.end(), mzInitCond.begin(), mzInitCond.end());

    // This zero is the (N+1)th spin on the RHS of the chain
    _mxStartVal.push_back(0);
    _myStartVal.push_back(0);
    _mzStartVal.push_back(0);

    /*
    int count2 = 0;
    for (int i =0; i < GV.GetNumSpins() + 2; i++) {
        if (++count2 % 10 == 0)
            std::cout << std::setw(12) << _mzStartVal[i] << std::endl;
        else
            std::cout << std::setw(12) << _mzStartVal[i] << ", ";
    }
    std::cout << "\n\n";
    */
}
void Numerical_Methods_Class::GilbertVectorsBothSides() {
    LinspaceClass GilbertDampingLHS;
    LinspaceClass GilbertDampingRHS;

    std::vector<double> gilbertChain(_undampedNumSpins, _gilbertConst);


    GilbertDampingLHS.set_values(_gilbertUpper, _gilbertLower, _numGilbert, true);
    GilbertDampingRHS.set_values(_gilbertLower, _gilbertUpper, _numGilbert, true);
    std::vector<double> tempGilbertLHS = GilbertDampingLHS.generate_array();
    std::vector<double> tempGilbertRHS = GilbertDampingRHS.generate_array();


    _gilbertVector.insert(_gilbertVector.end(), tempGilbertLHS.begin(), tempGilbertLHS.end());
    _gilbertVector.insert(_gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());
    _gilbertVector.insert(_gilbertVector.end(), tempGilbertRHS.begin(), tempGilbertRHS.end());

    _gilbertVector.push_back(0);

    /*
    int count = 0;
    for (int i = 0; i < GV.GetNumSpins() + 2; i++) {
        if (++count % 10 == 0)
            std::cout << std::setw(12) << _gilbertVector[i] << std::endl;
        else
            std::cout << std::setw(12) << _gilbertVector[i] << ", ";
    }
    */
}
void Numerical_Methods_Class::SetupVectorsGilbert() {
    LinspaceClass GilbertDamping;
    // Handles the Gilbert Damping vector.
    if (_numGilbert < 0)
        { // Guard clause.
            std::cout << "numGilbert is less than zero!";
            exit(0);
        }

    if (_lhsDrive)
    {
        GilbertDamping.set_values(_gilbertLower, _gilbertUpper, _numGilbert, true);
        std::vector<double> tempGilbert = GilbertDamping.generate_array();

        std::vector<double> gilbertChain(GV.GetNumSpins() - _numGilbert, _gilbertConst);
        gilbertChain.insert(gilbertChain.end(), tempGilbert.begin(), tempGilbert.end());
        _gilbertVector.insert(_gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());

        _gilbertVector.push_back(0);
    }
    else
    {
        GilbertDamping.set_values(_gilbertUpper, _gilbertLower, _numGilbert, true);
        std::vector<double> scaledGilbertChain = GilbertDamping.generate_array();

        std::vector<double> gilbertChain(GV.GetNumSpins() - _numGilbert, _gilbertConst);
        _gilbertVector.insert(_gilbertVector.end(), scaledGilbertChain.begin(), scaledGilbertChain.end());
        _gilbertVector.insert(_gilbertVector.end(), gilbertChain.begin(), gilbertChain.end());

        _gilbertVector.push_back(0);
    }

    /* test code to find vector contents
    int count = 0;
    for (int i = 0; i < _gilbertVector.size(); i++) {
        if (++count % 10 == 0) {std::cout << std::setw(12) << _gilbertVector[i] << "\n";}
        else {std::cout << std::setw(12) << _gilbertVector[i] << " ";}
    }

    exit(0);
    */


}

void Numerical_Methods_Class::RK2() {
    
    SpinChainEigenSolverClass printtest;

    // Notifies the user of what code they are running
    std::cout << "\nYou are running the RK2 Spinchains code." << std::endl;

    // Create files to save the data. All files will have (namefile) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath()+"rk2_mx_"+GV.GetFileNameBase()+".csv");
    std::ofstream myRK2File(GV.GetFilePath()+"rk2_my_"+GV.GetFileNameBase()+".csv");
    std::ofstream mzRK2File(GV.GetFilePath()+"rk2_mz_"+GV.GetFileNameBase()+".csv");

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or halfstep (h/2)*/
    for (long iterationIndex = static_cast<long>(_startIterVal); iterationIndex <= (long) _stopIterVal; iterationIndex++) {

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
            double mx1 = _mxStartVal[spin], mx1LHS = _mxStartVal[LHS_spin], mx1RHS = _mxStartVal[RHS_spin];
            // The my components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double my1 = _myStartVal[spin], my1LHS = _myStartVal[LHS_spin], my1RHS = _myStartVal[RHS_spin];
            // The mz components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mz1 = _mzStartVal[spin], mz1LHS = _mzStartVal[LHS_spin], mz1RHS = _mzStartVal[RHS_spin];

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

        // The estimate of the mx value for the next iteration of iterationIndex calculated using the RK2 Midpoint rule
        std::vector<double> mxNextVal(GV.GetNumSpins()+2,0);
        // The estimate of the my value for the next iteration of iterationIndex
        std::vector<double> myNextVal(GV.GetNumSpins()+2,0);
        // The estimate of the mz value for the next iteration of iterationIndex
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

            mxNextVal[spin] = _mxStartVal[spin] + mx2K2*_stepsize;
            myNextVal[spin] = _myStartVal[spin] + my2K2*_stepsize;
            mzNextVal[spin] = _mzStartVal[spin] + mz2K2*_stepsize;

        } // Final line of the RK2 solvers for this iteration. Everything below here is part of the class function, but not the internal RK2 stage loops

        // Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear these between
        // loop iterations sometimes led to incorrect values cropping up
        _mxStartVal.clear();
        _myStartVal.clear();
        _mzStartVal.clear();
        mxEstMid.clear();
        myEstMid.clear();
        mzEstMid.clear();

        /* Output function to write magnetic moment components to the terminal and/or files. Modulus component of IF
         * statement (default: 0.01 indicates how often the writing should occur. A value of 0.01 would mean writing
         * should occur every 1% of progress through the simulation*/
        if ( iterationIndex % int(_stopIterVal*0.01) == 0 ) {
            std::cout << "Reporting Point: " << iterationIndex << " iterations." << std::endl;
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
        _mxStartVal = mxNextVal;
        _myStartVal = myNextVal;
        _mzStartVal = mzNextVal;
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();
    myRK2File.close();
    mzRK2File.close();

    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "Finished RK2 with: stepSize = " << _stepsize << "; itermax = " << _stopIterVal << "; filename = " << GV.GetFileNameBase() <<  std::endl;
}

void Numerical_Methods_Class::RK2LLG() {

    progressbar bar(100);

    // Notifies the user of what code they are running.
    InformUserOfCodeType();

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath()+"rk2_mx_"+GV.GetFileNameBase()+".csv");
    //std::ofstream myRK2File(GV.GetFilePath()+"rk2_my_"+GV.GetFileNameBase()+".csv");
    //std::ofstream mzRK2File(GV.GetFilePath()+"rk2_mz_"+GV.GetFileNameBase()+".csv");

    CreateFileHeader(mxRK2File, _saveAllSpins, _onlyShowFinalState);
    //CreateFileHeader(myRK2File, _saveAllSpins, _onlyShowFinalState);
    //CreateFileHeader(mzRK2File, _saveAllSpins, _onlyShowFinalState);

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or half-step (h/2)*/
    for (int iterationIndex = _startIterVal; iterationIndex <= _stopIterVal; iterationIndex++) {
        // Doesn't work on Windows
        if (iterationIndex % (_stopIterVal / 100) == 0)
            {
                bar.update();
            }

        SetShockwaveConditions(iterationIndex);

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
            double mx1 = _mxStartVal[spin], mx1LHS = _mxStartVal[spinLHS], mx1RHS = _mxStartVal[spinRHS];
            // The my components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double my1 = _myStartVal[spin], my1LHS = _myStartVal[spinLHS], my1RHS = _myStartVal[spinRHS];
            // The mz components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mz1 = _mzStartVal[spin], mz1LHS = _mzStartVal[spinLHS], mz1RHS = _mzStartVal[spinRHS];

            double mxK1, myK1, mzK1; // These are the estimations of the slopes at the beginning of the interval for each magnetic moment component
            double heffXK1, heffYK1, heffZK1; // The effective field component acting upon each spin

            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                heffXK1 = _chainJVals[spinLHS] * mx1LHS + _chainJVals[spin] * mx1RHS +
                          _biasFieldDriving * cos(_drivingAngFreq * t0);
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

        // The estimate of the mx value for the next iteration of iterationIndex calculated using the RK2 Midpoint rule
        std::vector<double> mxNextVal(GV.GetNumSpins() + 2,0);
        // The estimate of the mY value for the next iteration of iterationIndex
        std::vector<double> myNextVal(GV.GetNumSpins() + 2,0);
        // The estimate of the mz value for the next iteration of iterationIndex
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
                heffXK2 = _chainJVals[spinLHS] * mx2LHS + _chainJVals[spin] * mx2RHS
                          + _biasFieldDriving * cos(_drivingAngFreq * t0HalfStep);
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

            mxNextVal[spin] = _mxStartVal[spin] + mxK2 * _stepsize;
            myNextVal[spin] = _myStartVal[spin] + myK2 * _stepsize;
            mzNextVal[spin] = _mzStartVal[spin] + mzK2 * _stepsize;
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        // Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
        // these between loop iterations sometimes led to incorrect values cropping up.
        _mxStartVal.clear();
        _myStartVal.clear();
        _mzStartVal.clear();
        mxEstMid.clear();
        myEstMid.clear();
        mzEstMid.clear();

        SaveDataToFile(_saveAllSpins, mxRK2File, mxNextVal, iterationIndex, _onlyShowFinalState);
        //SaveDataToFile(_saveAllSpins, myRK2File, myNextVal, iterationIndex, _onlyShowFinalState);
        //SaveDataToFile(_saveAllSpins, mzRK2File, mzNextVal, iterationIndex, _onlyShowFinalState);

        /* Sets the final value of the current iteration of the loop (y_(n+1) in textbook's notation) to be the starting
         * value of the next iteration (y_n) */
        _mxStartVal = mxNextVal;
        _myStartVal = myNextVal;
        _mzStartVal = mzNextVal;
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();
    //myRK2File.close();
    //mzRK2File.close();
    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::InformUserOfCodeType() {
    // Informs the user of the code type they are running, including: solver type; special modules.
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

    outputFileName << "Key Data\n" << std::endl;
    outputFileName << "Bias Field (H0) [T], Bias Field (Driving) [T], "
                          "Bias Field Driving Scale, Driving Frequency [Hz], Driving Region Start Site, Driving Region End Site, Driving Region Width,"
                          "Max. Sim. Time [s], Max. Exchange Val [T], Max. Iterations, Min. Exchange Val [T], "
                          "Num. DataPoints, Num. Spins, Stepsize (h), Damped Region Width\n";

    outputFileName << GV.GetBiasField() << ", " << _biasFieldDriving << ", " << _shockwaveScaling << ", "
                   << _drivingFreq << ", " << _drivingRegionLHS << ", " << _drivingRegionRHS << ", "
                   << _drivingRegionWidth << ", " << _maxSimTime << ", " << GV.GetExchangeMaxVal() << ", "
                   << _stopIterVal << ", " << GV.GetExchangeMinVal() << ", " << _numberOfDataPoints << ", "
                   << GV.GetNumSpins() << ", " << _stepsize << ", "<< _numGilbert << "\n\n";

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments );
    outputFileName << "Note(s):," << notesComments; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "\n[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n\n";

    CreateColumnHeaders(outputFileName, areAllSpinBeingSaved, onlyShowFinalState);

    std::cout << "\n";
}
void Numerical_Methods_Class::CreateColumnHeaders(std::ofstream &outputFileName, bool &areAllSpinBeingSaved, bool &onlyShowFinalState) {
    /* Creates the column headers for each spin site simulated. This code can change often, so compartmentalising it in
     * a separate function is necessary to reduce bugs.
     */
    if (areAllSpinBeingSaved or onlyShowFinalState) {
        // Print column heading for every spin simulated.
        outputFileName << "Time";
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
        if (current_iteration >= _stopIterVal * _iterToBeginShockwave)
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
        if (iteration % (_stopIterVal / _numberOfDataPoints) == 0) {
        //if (iteration == _stopIterVal) {
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
            if (iteration % (_stopIterVal / _numberOfDataPoints) == 0)
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

void Numerical_Methods_Class::StreamToString() {
    
    // Create an output string stream
    std::ostringstream stepsizeObj;
    std::ostringstream stopiterObj;
    
    stepsizeObj << _stepsize;
    stopiterObj << _stopIterVal;
    
    // Get string from output string stream
    _stepsizeString = stepsizeObj.str();
    _stopIterString = stopiterObj.str();
    
    stepsizeObj.clear();
    stopiterObj.clear();
}
void Numerical_Methods_Class::DebugOptions(std::vector<double> mxNextVal, std::vector<double> myNextVal, std::vector<double> mzNextVal, int spin, long iterationIndex) {
        if (mxNextVal[spin] >= 5000) {
            std::cout << "Error. Value of mx was greater than 5000 at spin(" << spin << "), iter("
                      << iterationIndex << ")." << std::flush;
            std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                      << _drivingRegionLHS << "; itermax = " << _stopIterVal << "; stepSize: "
                      << _stepsize << std::endl;
            exit(3);
        }

        if (myNextVal[spin] >= 5000) {
            std::cout << "Error. Value of my was greater than 5000 at spin(" << spin << "), iter("
                      << iterationIndex << ")." << std::flush;
            std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                      << _drivingRegionLHS << "; itermax = " << _stopIterVal << "; stepSize: "
                      << _stepsize << std::endl;
            exit(3);
        }

        if (mzNextVal[spin] >= 5000) {
            std::cout << "Error. Value of mz was greater than 5000 at spin(" << spin << "), iter("
                      << iterationIndex << ")." << std::flush;
            std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                      << _drivingRegionLHS << "; itermax = " << _stopIterVal << "; stepSize: "
                      << _stepsize << std::endl;
            exit(3);
        }
}
