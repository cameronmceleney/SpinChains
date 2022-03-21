#include "Numerical_Methods_Class.h"
#include "linspace.h"
#include "SpinChainEigenSolverClass.h"

void Numerical_Methods_Class::NMSetup() {

    _drivingFreq = 3.5 * 1e9;
    _drivingAngFreq = 2 * M_PI * _drivingFreq;
    _biasFieldDriving = 0.035;

    //std::cout << "Enter the LHS spin position for the driving region (minval=1): ";
    //std::cin >> _drivingRegionLHS;
    _drivingRegionLHS = 1;
    _drivingRegionWidth = int(GV.GetNumSpins() * 0.05);
    _drivingRegionRHS = _drivingRegionLHS + _drivingRegionWidth;


    // --- For testing on 20 Mar 22 ---
    //_drivingRegionLHS = 1;
    //_drivingRegionWidth = int(GV.GetNumSpins() * 0.05);
    //_drivingRegionRHS = 2;

    // _linearFMR = (_gyroMagConst / 2 * M_PI) * sqrt(_biasField * (_biasField + 4 * M_PI * _magSat)) / 1e9; // Only here for testing
    // --------------------------------

    _stepsize = 1e-15; // This should be at least (1 / _drivingFreq)
    _stepsizeHalf = _stepsize / 2.0;

    _stopIterVal = static_cast<int>(7e7);
    _numberOfDataPoints = _stopIterVal; // only here while I save every datapoint
    _maxSimTime = _stepsize * _stopIterVal;

    _hasShockwave = false;
    _iterToBeginShockwave = 0.5; // Value should be between [0.0, 1.0] inclusive.
    _shockwaveScaling = 5.0;
    _useLLG = true;
    _saveAllSpins = true;

    if (_drivingRegionRHS > GV.GetNumSpins()) {
        std::cout << "The width of the domain takes it past the maximum number of spins. Exiting...";
        //exit(3);
    }

    SetupVectors(); // Calls private class method to generate vectors needed for RK methods..
}
void Numerical_Methods_Class::SetupVectors() {

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

    // Notifies the user of what code they are running
    InformUserOfCodeType();

    // Create files to save the data. All files will have (namefile) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath()+"rk2_mx_"+GV.GetFileNameBase()+".csv");

    CreateFileHeader(mxRK2File, _saveAllSpins);

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or half-step (h/2)*/
    for (int iterationIndex = _startIterVal; iterationIndex <= _stopIterVal; iterationIndex++) {

        if (_hasShockwave and not _isShockwaveAlreadyOn) {
                if (iterationIndex >= _stopIterVal * _iterToBeginShockwave) {
                    // Shockwave begins once simulation is 50% complete
                    SetShockwaveConditions();
                    _isShockwaveAlreadyOn = true;
                }
        }

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
            } else if (spin > _drivingRegionRHS) {
                // The ELSE IF statement includes all spins along x which are not within the driving region
                heffXK1 = _chainJVals[spinLHS] * mx1LHS + _chainJVals[spin] * mx1RHS;
            } else {
                std::cout << "Error with selecting driving region spins in K1." << std::endl;
            }

            // No changes are made to the effective field in the y-direction
            heffYK1 = _chainJVals[spinLHS] * my1LHS + _chainJVals[spin] * my1RHS;
            // The bias field is applied in the z-direction and so it contributes to the effective field in the z-direction
            heffZK1 = _chainJVals[spinLHS] * mz1LHS + _chainJVals[spin] * mz1RHS + _biasField;

            if (_useLLG) {
                /* The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters
                 * for the first stage of RK2. */
                mxK1 = _gyroMagConst * (- (_gilbertConst * heffYK1 * mx1 * my1) + heffYK1 * mz1 - heffZK1 * (my1 + _gilbertConst*mx1*mz1) + _gilbertConst * heffXK1 * (pow(my1,2) + pow(mz1,2)));
                myK1 = _gyroMagConst * (-(heffXK1*mz1) + heffZK1 * (mx1 - _gilbertConst * my1 * mz1) + _gilbertConst * (heffYK1 * pow(mx1,2) - heffXK1 * mx1 * my1 + heffYK1 * pow(mz1,2)));
                mzK1 = _gyroMagConst * (heffXK1 * my1 + _gilbertConst * heffZK1*(pow(mx1,2) + pow(my1,2)) - _gilbertConst*heffXK1*mx1*mz1 - heffYK1 * (mx1 + _gilbertConst * my1 * mz1));
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
            else if (spin > _drivingRegionRHS)
                heffXK2 = _chainJVals[spinLHS] * mx2LHS + _chainJVals[spin] * mx2RHS;
            else
                std::cout << "Error with selecting driving region spins in K2." << std::endl;

            heffYK2 = _chainJVals[spinLHS] * my2LHS + _chainJVals[spin] * my2RHS;
            heffZK2 = _chainJVals[spinLHS] * mz2LHS + _chainJVals[spin] * mz2RHS + _biasField;

            if (_useLLG) {
                /* The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters
                 * for the second stage of RK2. */
                mxK2 = _gyroMagConst * (- (_gilbertConst * heffYK2 * mx2 * my2) + heffYK2 * mz2 - heffZK2 * (my2 + _gilbertConst*mx2*mz2) + _gilbertConst * heffXK2 * (pow(my2,2) + pow(mz2,2)));
                myK2 = _gyroMagConst * (-(heffXK2*mz2) + heffZK2 * (mx2 - _gilbertConst * my2 * mz2) + _gilbertConst * (heffYK2 * pow(mx2,2) - heffXK2 * mx2 * my2 + heffYK2 * pow(mz2,2)));
                mzK2 = _gyroMagConst * (heffXK2 * my2 + _gilbertConst * heffZK2*(pow(mx2,2) + pow(my2,2)) - _gilbertConst*heffXK2*mx2*mz2 - heffYK2 * (mx2 + _gilbertConst * my2 * mz2));
            } else {
                /* The magnetic moment components' coupled equations (obtained from the torque equation) with the
                 * parameters for the second stage of RK2. */
                mxK2 = -1 * _gyroMagConst * (my2 * heffZK2 - mz2 * heffYK2);
                myK2 = _gyroMagConst * (mx2 * heffZK2 - mz2 * heffXK2);
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

        SaveDataToFile(_saveAllSpins, mxRK2File, mxNextVal, iterationIndex);

        /* Sets the final value of the current iteration of the loop (y_(n+1) in textbook's notation) to be the starting
         * value of the next iteration (y_n) */
        _mxStartVal = mxNextVal;
        _myStartVal = myNextVal;
        _mzStartVal = mzNextVal;
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "\nFile can be found at:\n" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void Numerical_Methods_Class::RK2Shockwaves() {

    SpinChainEigenSolverClass printtest;

    // Sets the values of the driving field for before (Init) and after (Shock) the shockwave point respectively
    _biasFieldDrivingInit = _biasFieldDriving;
    _biasFieldDrivingShock = _biasFieldDriving * _biasFieldDrivingScale;
    _hasShockWaveBegan = false;

    // Notifies the user of what code they are running
    std::cout << "\nYou are running the RK2 Shockwave Spinchains code." << std::endl;
    // Create files to save the data. All files will have (FileNameBase) in them to make them clearly identifiable.
    std::ofstream mxRK2ShockwaveFile(GV.GetFilePath()+"rk2Shockwave_"+GV.GetFileNameBase()+".csv");

    CreateFileHeader(mxRK2ShockwaveFile, _saveAllSpins);

    std::cout << "\nBeginning simulation...";
    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or half-step (h/2)*/
    for (long iterationIndex = static_cast<long>(_startIterVal); iterationIndex <= (long) _stopIterVal; iterationIndex++) {


        if (iterationIndex >= (long)_stopIterVal*0.5){
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
        for (int spinMid1 = 1; spinMid1 <= GV.GetNumSpins()+1; spinMid1++) {
            /* The first stage is based upon finding the value of the slope at the beginning of the interval (k1). This
             * stage takes the start conditions as an input, and substitutes them into the LLG equation. */
            int spinToLHS1 = spinMid1 - 1, spinToRHS1 = spinMid1 + 1;
            //std::cout << "Spin #" << spinMid1 << " is using biasField " << _biasFieldDrivingUse << std::endl;
            // The mx components for the first stage for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mx1 = _mxStartVal[spinMid1], mx1LHS = _mxStartVal[spinToLHS1], mx1RHS = _mxStartVal[spinToRHS1];
            // The my components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double my1 = _myStartVal[spinMid1], my1LHS = _myStartVal[spinToLHS1], my1RHS = _myStartVal[spinToRHS1];
            // The mz components for the first  for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mz1 = _mzStartVal[spinMid1], mz1LHS = _mzStartVal[spinToLHS1], mz1RHS = _mzStartVal[spinToRHS1];

            double mx1K1, my1K1, mz1K1; // These are the estimations of the slopes at the beginning of the interval for each magnetic moment component
            double HeffX1K1, HeffY1K1, HeffZ1K1; // The effective field component acting upon each spin

            if (spinMid1 >= _drivingRegionLHS && spinMid1 <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                HeffX1K1 = _chainJVals[spinToLHS1] * mx1LHS + _chainJVals[spinMid1] * mx1RHS + _biasFieldDrivingUse*cos(_drivingAngFreq * t0);
            } else if (spinMid1 > _drivingRegionRHS) {
                // The ELSE IF statement includes all spins along x which are not within the driving region
                HeffX1K1 = _chainJVals[spinToLHS1] * mx1LHS + _chainJVals[spinMid1] * mx1RHS;
            } else {
                // No statement
            }
            // No changes are made to the effective field in the y-direction
            HeffY1K1 = _chainJVals[spinToLHS1] * my1LHS + _chainJVals[spinMid1] * my1RHS;
            // The bias field is applied in the z-direction and so it contributes to the effective field in the z-direction
            HeffZ1K1 = _chainJVals[spinToLHS1] * mz1LHS + _chainJVals[spinMid1] * mz1RHS + _biasField;

            /* The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the
             * first stage of RK2.*/
            mx1K1 = -1 * _gyroMagConst * (my1 * HeffZ1K1 - mz1 * HeffY1K1);
            my1K1 = +1 * _gyroMagConst * (mx1 * HeffZ1K1 - mz1 * HeffX1K1);
            mz1K1 = -1 * _gyroMagConst * (mx1 * HeffY1K1 - my1 * HeffX1K1);

            mxEstMid[spinMid1] = mx1 + mx1K1*_stepsizeHalf;
            myEstMid[spinMid1] = my1 + my1K1*_stepsizeHalf;
            mzEstMid[spinMid1] = mz1 + mz1K1*_stepsizeHalf;
        }

        // The estimate of the mx value for the next iteration of iterationIndex calculated using the RK2 Midpoint rule
        std::vector<double> mxNextVal(GV.GetNumSpins()+2, 0);
        // The estimate of the my value for the next iteration of iterationIndex
        std::vector<double> myNextVal(GV.GetNumSpins()+2, 0);
        // The estimate of the mz value for the next iteration of iterationIndex
        std::vector<double> mzNextVal(GV.GetNumSpins()+2, 0);

        for (int spinMid2 = 1; spinMid2 <= GV.GetNumSpins()+1; spinMid2++) {
            /* The second stage uses the previously found k1 value, as well as the initial conditions, to determine the
             * value of the slope (k2) at the midpoint. In RK2, the values of k1 and k2 can then be jointly used to
             * estimate the next point of the function through a weighted average of k1 & k2.
             *
             * In this loop the definitions of variables follow a similar format to Stage1.*/

            int spinToLHS2 = spinMid2 - 1, spinToRHS2 = spinMid2 + 1;
            double mx2 = mxEstMid[spinMid2], mx2LHS = mxEstMid[spinToLHS2], mx2RHS = mxEstMid[spinToRHS2];
            double my2 = myEstMid[spinMid2], my2LHS = myEstMid[spinToLHS2], my2RHS = myEstMid[spinToRHS2];
            double mz2 = mzEstMid[spinMid2], mz2LHS = mzEstMid[spinToLHS2], mz2RHS = mzEstMid[spinToRHS2];

            double mx2K2, my2K2, mz2K2;
            double HeffX2K2, HeffY2K2, HeffZ2K2;

            if (spinMid2 >= _drivingRegionLHS && spinMid2 <= _drivingRegionRHS) {
                // If a spin is driven during Stage 1 of an RK method, then it must be driven throughout the rest of the method's stages
                HeffX2K2 = _chainJVals[spinToLHS2] * mx2LHS + _chainJVals[spinMid2] * mx2RHS + _biasFieldDrivingUse*cos(_drivingAngFreq * t0HalfStep);
            } else if (spinMid2 > _drivingRegionRHS) {
                HeffX2K2 = _chainJVals[spinToLHS2] * mx2LHS + _chainJVals[spinMid2] * mx2RHS;
            } else {
                // No statement
            }

            HeffY2K2 = _chainJVals[spinToLHS2] * my2LHS + _chainJVals[spinMid2] * my2RHS;
            HeffZ2K2 = _chainJVals[spinToLHS2] * mz2LHS + _chainJVals[spinMid2] * mz2RHS + _biasField;

            mx2K2 = -1 * _gyroMagConst * ( my2 * HeffZ2K2 - mz2 * HeffY2K2 );
            my2K2 = +1 * _gyroMagConst * ( mx2 * HeffZ2K2 - mz2 * HeffX2K2 );
            mz2K2 = -1 * _gyroMagConst * ( mx2 * HeffY2K2 - my2 * HeffX2K2 );

            mxNextVal[spinMid2] = _mxStartVal[spinMid2] + mx2K2*_stepsize;
            myNextVal[spinMid2] = _myStartVal[spinMid2] + my2K2*_stepsize;
            mzNextVal[spinMid2] = _mzStartVal[spinMid2] + mz2K2*_stepsize;
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
         * should occur every 1% of progress through the simulation)*/
        if ( iterationIndex % int(_stopIterVal*(1.0/_numberOfDataPoints)) == 0 ) { // Value MUST be 1.0 to ensure correct casting

            //
            mxRK2ShockwaveFile << mxNextVal[_drivingRegionLHS] << ", "
                               << mxNextVal[_drivingRegionRHS] << ", "
                               << mxNextVal[int(1+GV.GetNumSpins()*0.5)]
                               << ", " << mxNextVal[int(1+GV.GetNumSpins())] << "\n";

        }

        /* Sets the final value of the current iteration of the loop (y_(n+1) in textbook's notation) to be the starting
         * value of the next iteration (y_n) */
        _mxStartVal = mxNextVal;
        _myStartVal = myNextVal;
        _mzStartVal = mzNextVal;

        mxNextVal.clear();
        myNextVal.clear();
        mzNextVal.clear();
    } // Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2ShockwaveFile.close();

    // Provides key parameters to user for their log. Filename can be copy/pasted from terminal to a plotter function in Python
    std::cout << "\nFinished RK2 with: stepsize = " << _stepsize << "; itermax = " << _stopIterVal << "; filename = " << GV.GetFileNameBase() <<  std::endl;
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
void Numerical_Methods_Class::CreateFileHeader(std::ofstream &outputFileName, bool &areAllSpinBeingSaved) {

    outputFileName << "Key Data\n" << std::endl;
    outputFileName << "Bias Field (H0) [T], Bias Field (Driving) [T], "
                          "Bias Field Driving Scale, Driving Frequency [Hz], Driving Region Start Site, Driving Region End Site, Driving Region Width,"
                          "Max. Sim. Time [s], Max. Exchange Val [T], Max. Iterations, Min. Exchange Val [T], "
                          "Num. DataPoints, Num. Spins, Stepsize (h)\n";

    outputFileName << _biasField << ", " << _biasFieldDriving << ", " << _biasFieldDrivingScale << ", "
                   << _drivingFreq << ", " << _drivingRegionLHS << ", " << _drivingRegionRHS - 1 << ", "
                   << _drivingRegionWidth << ", " << _maxSimTime << ", " << GV.GetExchangeMaxVal() << ", "
                   << _stopIterVal << ", " << GV.GetExchangeMinVal() << ", " << _numberOfDataPoints << ", "
                   << GV.GetNumSpins() << ", " << _stepsize << "\n\n";

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments );
    outputFileName << "Note(s):," << notesComments; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "\n[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n\n";
    // outputFileName << _drivingRegionLHS << ", " << _drivingRegionRHS - 1 << ", " << (GV.GetNumSpins()/2) << ", " << GV.GetNumSpins() << std::endl;
    if (areAllSpinBeingSaved) {
        // Print column heading for every spin simulated.
        outputFileName << "Time";
        for (int i = 1; i <= GV.GetNumSpins(); i++) {
            outputFileName << "," << i;
        }
        outputFileName << std::endl;

    } else {
        outputFileName << _drivingRegionLHS << ", " << _drivingRegionRHS - 1 << ", " << 3 << ", " << GV.GetNumSpins() << std::endl;
    }

}
void Numerical_Methods_Class::SetShockwaveConditions() {

    // If function is triggered, then the applied biasFieldDriving is increased by the scale factor _shockwaveScaling
    _biasFieldDriving *= _shockwaveScaling;

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

void Numerical_Methods_Class::SaveDataToFile(bool &_areAllSpinBeingSaved, std::ofstream &outputFileName,
                                             std::vector<double> &arrayToWrite, int iteration) {

    if (_areAllSpinBeingSaved) {
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
    } else {
            outputFileName << arrayToWrite[_drivingRegionLHS] << ", "
                           << arrayToWrite[_drivingRegionRHS] << ", "
                           << arrayToWrite[static_cast<int>(3 + 1)] << ", "
                           << arrayToWrite[static_cast<int>(1 + GV.GetNumSpins())] << "\n";

    }

}
