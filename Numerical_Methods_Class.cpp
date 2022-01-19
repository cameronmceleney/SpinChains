#include "Numerical_Methods_Class.h"
#include "linspace.h"
#include "SpinChainEigenSolverClass.h"

void Numerical_Methods_Class::NMSetup() {

    _drivingAngularFrequency = 2 * M_PI * _drivingFrequency;

    std::cout << "Enter the LHS spin position for the driving region: ";
    std::cin >> _drivingRegionLHS;

    std::cout << "Enter the stepsize: ";
    std::cin >> _stepSize;
    _stepsizeHalf = _stepSize / 2;

    _drivingRegionWidth = GV.GetNumSpins() * 0.05;
    _drivingRegionRHS = _drivingRegionLHS + _drivingRegionWidth;

    std::cout << "Enter the maximum number of iterations: ";
    std::cin >> _stopIterationValue; // Can be inputted in scientific notation or as a float

    _maxSimulatedTime = _stepSize * _stopIterationValue;
    _initialMagMomentZ = _magnetisationSaturation;

    std::cout << "\nThis will simulate a time of " << _maxSimulatedTime << "[s]." << std::endl;
}

void Numerical_Methods_Class::RK2() {

    LinspaceClass SpinChainExchange;
    SpinChainEigenSolverClass printtest;

    // Notifies the user of what code they are running
    std::cout << "\nYou are running the RK2 tester chainspin code." << std::endl;
    //std::cout << "DrivingLHS: " << _drivingRegionLHS << "; drivingRHS: " << _drivingRegionRHS << "; driving width: " << _drivingRegionWidth << std::endl;
    if(_drivingRegionRHS > GV.GetNumSpins()) {
        std::cout << "The width of the domain takes it past the maximum number of spins. Exiting...";
        exit(3);
    }
    int c_numberOfSpinPairs = GV.GetNumSpins() - 1; // Used to indicate array lengths and tidy notation
    SpinChainExchange.set_values(GV.GetExchangeMinVal(), GV.GetExchangeMaxVal(), c_numberOfSpinPairs,true);
    SpinChainExchange.generate_array();
    _chainExchangeValues = SpinChainExchange.build_spinchain();

    //TODO Turn the initial conditions lines into a separate function
    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mXInitCond(GV.GetNumSpins(), _initialMagMomentX), mYInitCond(GV.GetNumSpins(), _initialMagMomentY), mZInitCond(GV.GetNumSpins(), _initialMagMomentZ);
    std::vector<double> mXEstStart{0}, mYEstStart{0}, mZEstStart{0}; // Magnetic Component (m), Axis (X), Estimate (Est), Initial Time (Start).
    // Appends initial conditions to the vectors
    mXEstStart.insert(mXEstStart.end(), mXInitCond.begin(), mXInitCond.end());
    mYEstStart.insert(mYEstStart.end(), mYInitCond.begin(), mYInitCond.end());
    mZEstStart.insert(mZEstStart.end(), mZInitCond.begin(), mZInitCond.end());
    // Delete temporary vectors that held the initial conditions
    mXInitCond.clear();
    mYInitCond.clear();
    mZInitCond.clear();

    // This zero is the (N+1)th spin on the RHS of the chain
    mXEstStart.push_back(0);
    mYEstStart.push_back(0);
    mZEstStart.push_back(0);

    // Create files to save the data. All files will have (namefile) in them to make them clearly identifiable.
    std::ofstream mXFile(GV.GetFilePath()+"rk2_mx_"+GV.GetFileNameBase()+".csv");
    std::ofstream mYFile(GV.GetFilePath()+"rk2_my_"+GV.GetFileNameBase()+".csv");
    std::ofstream mZFile(GV.GetFilePath()+"rk2_mz_"+GV.GetFileNameBase()+".csv");

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or halfstep (h/2)*/
    for (long iterationIndex = _startIterationValue; iterationIndex <= (long) _stopIterationValue; iterationIndex++) {

        _totalTime += _stepSize;
        double t0 = _totalTime; // The initial time of the iteration, and the time at the first stage; the start of the interval and the first step of RK2. Often called 't0' in literature
        double t0HalfStep = _totalTime + _stepsizeHalf; // The time at the midpoint of the interval; where the second stage of RK2 is. Often called 't0*h/2' in literature

        // The estimate of the slope for the x-axis magnetic moment component at the midpoint
        std::vector<double> mXEstMid(GV.GetNumSpins()+2, 0);
        // The estimate of the slope for the y-axis magnetic moment component at the midpoint
        std::vector<double> mYEstMid(GV.GetNumSpins()+2, 0);
        // The estimate of the slope for the z-axis magnetic moment component at the midpoint
        std::vector<double> mZEstMid(GV.GetNumSpins()+2, 0);

        // Loop the 0th and final spins as they will always be zero-valued
        for (int spin = 1; spin <= GV.GetNumSpins()+1; spin++) {
            /* The first stage is based upon finding the value of the slope at the beginning of the interval (k1). This
             * stage takes the start conditions as an input, and substitutes them into the LLG equation. */
            int LHS_spin = spin - 1, RHS_spin = spin + 1;

            // The mX components for the first stage for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mXStage1 = mXEstStart[spin], mXStage1LHS = mXEstStart[LHS_spin], mXStage1RHS = mXEstStart[RHS_spin];
            // The mY components for the first stage for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mYStage1 = mYEstStart[spin], mYStage1LHS = mYEstStart[LHS_spin], mYStage1RHS = mYEstStart[RHS_spin];
            // The mZ components for the first stage for the: current spin site; site to the left (LHS); site to the right (RHS)
            double mZStage1 = mZEstStart[spin], mZStage1LHS = mZEstStart[LHS_spin], mZStage1RHS = mZEstStart[RHS_spin];

            double mXStage1K1, mYStage1K1, mZStage1K1; // These are the estimations of the slopes at the beginning of the interval for each magnetic moment component
            double HeffXStage1K1, HeffYStage1K1, HeffZStage1K1; // The effective field component acting upon each spin

            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
                HeffXStage1K1 = _chainExchangeValues[LHS_spin] * mXStage1LHS + _chainExchangeValues[spin] * mXStage1RHS + _biasFieldDrivingAmplitude*cos(_drivingAngularFrequency * t0);

            } else {
                // The else statement includes all spins along x which are not within the driving region
                HeffXStage1K1 = _chainExchangeValues[LHS_spin] * mXStage1LHS + _chainExchangeValues[spin] * mXStage1RHS;
            }
            // No changes are made to the effective field in the y-direction
            HeffYStage1K1 = _chainExchangeValues[LHS_spin] * mYStage1LHS + _chainExchangeValues[spin] * mYStage1RHS;
            // The bias field is applied in the z-direction and so it contributes to the effective field in the z-direction
            HeffZStage1K1 = _chainExchangeValues[LHS_spin] * mZStage1LHS + _chainExchangeValues[spin] * mZStage1RHS + _biasField;

            /* The magnetic moment components' coupled equations (obtained from LLG equation) with the parameters for the
             * first stage of RK2.*/
            mXStage1K1 = -1 * _gyroscopicMagneticConstant * (mYStage1 * HeffZStage1K1 - mZStage1 * HeffYStage1K1);
            mYStage1K1 = +1 * _gyroscopicMagneticConstant * (mXStage1 * HeffZStage1K1 - mZStage1 * HeffXStage1K1);
            mZStage1K1 = -1 * _gyroscopicMagneticConstant * (mXStage1 * HeffYStage1K1 - mYStage1 * HeffXStage1K1);

            mXEstMid[spin] = mXStage1 + mXStage1K1*_stepsizeHalf;
            mYEstMid[spin] = mYStage1 + mYStage1K1*_stepsizeHalf;
            mZEstMid[spin] = mZStage1 + mZStage1K1*_stepsizeHalf;
        }

        // The estimate of the mX value for the next iteration of iterationIndex calculated using the RK2 Midpoint rule
        std::vector<double> mXNextVal(GV.GetNumSpins()+2,0);
        // The estimate of the mY value for the next iteration of iterationIndex
        std::vector<double> mYNextVal(GV.GetNumSpins()+2,0);
        // The estimate of the mZ value for the next iteration of iterationIndex
        std::vector<double> mZNextVal(GV.GetNumSpins()+2,0);

        for (int spin = 1; spin <= GV.GetNumSpins()+1; spin++) {
            /* The second stage uses the previously found k1 value, as well as the initial conditions, to determine the
             * value of the slope (k2) at the midpoint. In RK2, the values of k1 and k2 can then be jointly used to
             * estimate the next point of the function through a weighted average of k1 & k2.
             * 
             * In this loop the definitions of variables follow a similar format to Stage1.*/

            int LHS_spin = spin - 1, RHS_spin = spin + 1;
            double mXStage2 = mXEstMid[spin], mXStage2LHS = mXEstMid[LHS_spin], mXStage2RHS = mXEstMid[RHS_spin];
            double mYStage2 = mYEstMid[spin], mYStage2LHS = mYEstMid[LHS_spin], mYStage2RHS = mYEstMid[RHS_spin];
            double mZStage2 = mZEstMid[spin], mZStage2LHS = mZEstMid[LHS_spin], mZStage2RHS = mZEstMid[RHS_spin];

            double mXStage2K2, mYStage2K2, mZStage2K2;
            double HeffXStage2K2, HeffYStage2K2, HeffZStage2K2;

            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                // Driving region must be consistently applied at every stage of the RK2 method
                HeffXStage2K2 = _chainExchangeValues[LHS_spin] * mXStage2LHS + _chainExchangeValues[spin] * mXStage2RHS + _biasFieldDrivingAmplitude*cos(_drivingAngularFrequency * t0HalfStep);
            } else {
                HeffXStage2K2 = _chainExchangeValues[LHS_spin] * mXStage2LHS + _chainExchangeValues[spin] * mXStage2RHS;
            }
            HeffYStage2K2 = _chainExchangeValues[LHS_spin] * mYStage2LHS + _chainExchangeValues[spin] * mYStage2RHS;
            HeffZStage2K2 = _chainExchangeValues[LHS_spin] * mZStage2LHS + _chainExchangeValues[spin] * mZStage2RHS + _biasField;

            mXStage2K2 = -1 * _gyroscopicMagneticConstant * ( mYStage2*HeffZStage2K2 - mZStage2*HeffYStage2K2 );
            mYStage2K2 = +1 * _gyroscopicMagneticConstant * ( mXStage2*HeffZStage2K2 - mZStage2*HeffXStage2K2 );
            mZStage2K2 = -1 * _gyroscopicMagneticConstant * ( mXStage2*HeffYStage2K2 - mYStage2*HeffXStage2K2 );

            mXNextVal[spin] = mXEstStart[spin] + mXStage2K2*_stepSize;
            mYNextVal[spin] = mYEstStart[spin] + mYStage2K2*_stepSize;
            mZNextVal[spin] = mZEstStart[spin] + mZStage2K2*_stepSize;


            if (_shouldDebug) {
                if (mXNextVal[spin] >= 5000) {
                    std::cout << "Error. Value of mx was greater than 5000 at spin(" << spin << "), iter("
                              << iterationIndex << ")." << std::flush;
                    std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                              << _drivingRegionLHS << "; itermax = " << _stopIterationValue << "; stepSize: "
                              << _stepSize << std::endl;
                    exit(3);
                }

                if (mYNextVal[spin] >= 5000) {
                    std::cout << "Error. Value of my was greater than 5000 at spin(" << spin << "), iter("
                              << iterationIndex << ")." << std::flush;
                    std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                              << _drivingRegionLHS << "; itermax = " << _stopIterationValue << "; stepSize: "
                              << _stepSize << std::endl;
                    exit(3);
                }

                if (mZNextVal[spin] >= 5000) {
                    std::cout << "Error. Value of mz was greater than 5000 at spin(" << spin << "), iter("
                              << iterationIndex << ")." << std::flush;
                    std::cout << " Test info as follows: numSpins = " << GV.GetNumSpins() << "; starting spin = "
                              << _drivingRegionLHS << "; itermax = " << _stopIterationValue << "; stepSize: "
                              << _stepSize << std::endl;
                    exit(3);
                }
            }

        }

        // Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp
        mXEstStart.clear();
        mYEstStart.clear();
        mZEstStart.clear();
        mXEstMid.clear();
        mYEstMid.clear();
        mZEstMid.clear();

        /* Output function to write magnetic moment components to the terminal and/or files. Modulus component of IF
        */
        if ( iterationIndex % int(_stopIterationValue*0.01) == 0 ) {
            std::cout << "Reporting Point: " << iterationIndex << " iterations." << std::endl;
            /*std::cout << "Numspins: " << GV.GetNumSpins() << "; const Jmin: " << GV.GetExchangeMinVal() << "; const Jmax: " << GV.GetExchangeMaxVal() << std::endl;
            std::cout << "RegionLHS: " << _drivingRegionLHS << "; RegionWidth: " << _drivingRegionWidth << "; RegionRHS: " << _drivingRegionRHS << std::endl;
            std::cout << "StepSize: " << _stepSize << "; HalfStepSize: " << _stepsizeHalf << "; TotalTime: " << _totalTime << "\n\n" << std::endl;
            */
            //printtest.PrintVector(mXNextVal);
            for (int j = 1; j < GV.GetNumSpins() + 1; j++) {
                mXFile << mXNextVal[j] << ",";
                mYFile << mYNextVal[j] << ",";
                mZFile << mZNextVal[j] << ",";
                if (j == GV.GetNumSpins()) {
                    mXFile << mXNextVal[j] << std::flush;
                    mYFile << mYNextVal[j] << std::flush;
                    mZFile << mZNextVal[j] << std::flush;
                }
            }
            mXFile << std::endl;
            mYFile << std::endl;
            mZFile << std::endl;
        }

        mXEstStart = mXNextVal;
        mYEstStart = mYNextVal;
        mZEstStart = mZNextVal;

    }
    mXFile.close();
    mYFile.close();
    mZFile.close();

    std::cout << "Finished with: stepSize = " << _stepSize << "; itermax = " << _stopIterationValue << "; filename = " << GV.GetFileNameBase() <<  std::endl;

}

void Numerical_Methods_Class::StreamToString() {
    // Create an output string stream
    std::ostringstream stepsizeObj;
    std::ostringstream stopiterObj;
    stepsizeObj << _stepSize;
    stopiterObj << _stopIterationValue;
    // Get string from output string stream
    _stepsizeString = stepsizeObj.str();
    _stopIterString = stopiterObj.str();
    stepsizeObj.clear();
    stopiterObj.clear();
}