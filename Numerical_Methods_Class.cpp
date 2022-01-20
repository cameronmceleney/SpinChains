#include "Numerical_Methods_Class.h"
#include "linspace.h"
#include "SpinChainEigenSolverClass.h"

void Numerical_Methods_Class::NMSetup() {

    LinspaceClass SpinChainExchange;
    
    // Calculate Numerical_Methods_Class members which are only dependent upon the user's prior input
    _drivingAngFreq = 2 * M_PI * _drivingFreq;
    _mzInit = _magSat;
    _numberOfSpinPairs = GV.GetNumSpins() - 1;

    std::cout << "Enter the LHS spin position for the driving region (minval=1): ";
    std::cin >> _drivingRegionLHS;
    _drivingRegionWidth = GV.GetNumSpins() * 0.05;
    _drivingRegionRHS = _drivingRegionLHS + _drivingRegionWidth;

    std::cout << "Enter the stepsize: ";
    std::cin >> _stepsize;
    _stepsizeHalf = _stepsize / 2;
    
    std::cout << "Enter the maximum number of iterations: ";
    std::cin >> _stopIterVal; // Can be inputted in scientific notation or as a float
    _maxSimTime = _stepsize * _stopIterVal;

    std::cout << "\nThis will simulate a time of " << _maxSimTime << "[s]." << std::endl;

    if(_drivingRegionRHS > GV.GetNumSpins()) {
        std::cout << "The width of the domain takes it past the maximum number of spins. Exiting...";
        exit(3);
    }
    
    // The linearly spaced vector is saved as the class member '_chainJVals' simply to increase code readability
    SpinChainExchange.set_values(GV.GetExchangeMinVal(), GV.GetExchangeMaxVal(), _numberOfSpinPairs, true);
    SpinChainExchange.generate_array();
    _chainJVals = SpinChainExchange.build_spinchain();
    
    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mxInitCond(GV.GetNumSpins(), _mxInit), myInitCond(GV.GetNumSpins(), _myInit), mzInitCond(GV.GetNumSpins(), _mzInit);
    
    // Appends initial conditions to the vectors
    _mxStartVal.insert(_mxStartVal.end(), mxInitCond.begin(), mxInitCond.end());
    _myStartVal.insert(_myStartVal.end(), myInitCond.begin(), myInitCond.end());
    _mzStartVal.insert(_mzStartVal.end(), mzInitCond.begin(), mzInitCond.end());
    // Delete temporary vectors that held the initial conditions
    mxInitCond.clear();
    myInitCond.clear();
    mzInitCond.clear();
    // This zero is the (N+1)th spin on the RHS of the chain
    _mxStartVal.push_back(0);
    _myStartVal.push_back(0);
    _mzStartVal.push_back(0);
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
    for (long iterationIndex = _startIterVal; iterationIndex <= (long) _stopIterVal; iterationIndex++) {

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
            HeffZ1K1 = _chainJVals[LHS_spin] * mz1LHS + _chainJVals[spin] * mz1RHS + _biasField;

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
            HeffZ2K2 = _chainJVals[LHS_spin] * mz2LHS + _chainJVals[spin] * mz2RHS + _biasField;

            mx2K2 = -1 * _gyroMagConst * ( my2*HeffZ2K2 - mz2*HeffY2K2 );
            my2K2 = +1 * _gyroMagConst * ( mx2*HeffZ2K2 - mz2*HeffX2K2 );
            mz2K2 = -1 * _gyroMagConst * ( mx2*HeffY2K2 - my2*HeffX2K2 );

            mxNextVal[spin] = _mxStartVal[spin] + mx2K2*_stepsize;
            myNextVal[spin] = _myStartVal[spin] + my2K2*_stepsize;
            mzNextVal[spin] = _mzStartVal[spin] + mz2K2*_stepsize;

            if (_shouldDebug) {
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

void Numerical_Methods_Class::RK2Shockwaves() {

    SpinChainEigenSolverClass printtest;

    // Sets the values of the driving field for before (Init) and after (Shock) the shockwave point respectively
    _biasFieldDrivingInit = _biasFieldDriving;
    _biasFieldDrivingShock = _biasFieldDriving;
    _hasShockWaveBegan = false;
    // Notifies the user of what code they are running
    std::cout << "\nYou are running the RK2 Shockwave Spinchains code." << std::endl;

    // Create files to save the data. All files will have (FileNameBase) in them to make them clearly identifiable.
    std::ofstream mxRK2ShockwaveLHSDriving(GV.GetFilePath()+"rk2Shockwave_mxinputLHSdriving_"+GV.GetFileNameBase()+".csv");
    std::ofstream mxRK2ShockwaveRHSDriving(GV.GetFilePath()+"rk2Shockwave_mxinputRHSdriving_"+GV.GetFileNameBase()+".csv");
    std::ofstream mxRK2ShockwaveMid(GV.GetFilePath()+"rk2Shockwave_mxoutputmidpoint_"+GV.GetFileNameBase()+".csv");
    std::ofstream mxRK2ShockwaveEnd(GV.GetFilePath()+"rk2Shockwave_mxoutputendpoint_"+GV.GetFileNameBase()+".csv");

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or halfstep (h/2)*/
    for (long iterationIndex = _startIterVal; iterationIndex <= (long) _stopIterVal; iterationIndex++) {

        if (iterationIndex >= _stopIterVal*0.5){
            // Shockwave begins once simulation is 50% complete
            _hasShockWaveBegan = true;
        }

        if (_hasShockWaveBegan){
            // Changes used driving field to be post-shockwave value
            _biasFieldDrivingUse = _biasFieldDrivingShock;
        } else if (!_hasShockWaveBegan) {
            // Keeps used driving field value to be pre-shockwave value
            _biasFieldDrivingUse = _biasFieldDrivingInit;
        } else {
            // No statement
        }

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
                HeffX1K1 = _chainJVals[LHS_spin] * mx1LHS + _chainJVals[spin] * mx1RHS + _biasFieldDrivingUse*cos(_drivingAngFreq * t0);
            } else {
                // The else statement includes all spins along x which are not within the driving region
                HeffX1K1 = _chainJVals[LHS_spin] * mx1LHS + _chainJVals[spin] * mx1RHS;
            }

            // No changes are made to the effective field in the y-direction
            HeffY1K1 = _chainJVals[LHS_spin] * my1LHS + _chainJVals[spin] * my1RHS;
            // The bias field is applied in the z-direction and so it contributes to the effective field in the z-direction
            HeffZ1K1 = _chainJVals[LHS_spin] * mz1LHS + _chainJVals[spin] * mz1RHS + _biasField;

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
                HeffX2K2 = _chainJVals[LHS_spin] * mx2LHS + _chainJVals[spin] * mx2RHS + _biasFieldDrivingUse*cos(_drivingAngFreq * t0HalfStep);
            } else {
                HeffX2K2 = _chainJVals[LHS_spin] * mx2LHS + _chainJVals[spin] * mx2RHS;
            }

            HeffY2K2 = _chainJVals[LHS_spin] * my2LHS + _chainJVals[spin] * my2RHS;
            HeffZ2K2 = _chainJVals[LHS_spin] * mz2LHS + _chainJVals[spin] * mz2RHS + _biasField;

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
        if ( iterationIndex % int(_stopIterVal*0.001) == 0 ) {
            std::cout << "Reporting Point: " << iterationIndex << " iterations. Shockwave bool: " << _hasShockWaveBegan << " with driving val: " << _biasFieldDrivingUse <<std::endl;
            /*  Code this is useful for debugging
             *  std::cout << "Numspins: " << GV.GetNumSpins() << "; const Jmin: " << GV.GetExchangeMinVal() << "; const Jmax: " << GV.GetExchangeMaxVal() << std::endl;
             *  std::cout << "RegionLHS: " << _drivingRegionLHS << "; RegionWidth: " << _drivingRegionWidth << "; RegionRHS: " << _drivingRegionRHS << std::endl;
             *  std::cout << "StepSize: " << _stepsize << "; HalfStepSize: " << _stepsizeHalf << "; TotalTime: " << _totalTime << "\n\n" << std::endl;
             *  printtest.PrintVector(mxNextVal); */

            // This is the input signal. Saves all mx values at spin site 1 (spin site zero isn't really a spin).
            mxRK2ShockwaveLHSDriving << mxNextVal[_drivingRegionLHS] << std::endl;
            mxRK2ShockwaveRHSDriving << mxNextVal[_drivingRegionRHS] << std::endl;
            // This is the output signal. Saves all mx values at second last array element, and the final real spin site.
            mxRK2ShockwaveMid << mxNextVal[1+(GV.GetNumSpins()/2)] << std::endl;
            mxRK2ShockwaveEnd << mxNextVal[GV.GetNumSpins()-5] << std::endl;

            // Housekeeping
            //mxRK2ShockwaveStart << std::endl;
            //mxRK2ShockwaveEnd << std::endl;
        }

        /* Sets the final value of the current iteration of the loop (y_(n+1) in textbook's notation) to be the starting
         * value of the next iteration (y_n) */
        _mxStartVal = mxNextVal;
        _myStartVal = myNextVal;
        _mzStartVal = mzNextVal;
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2ShockwaveLHSDriving.close();
    mxRK2ShockwaveRHSDriving.close();
    mxRK2ShockwaveMid.close();
    mxRK2ShockwaveEnd.close();

    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "Finished RK2 with: stepSize = " << _stepsize << "; itermax = " << _stopIterVal << "; filename = " << GV.GetFileNameBase() <<  std::endl;
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