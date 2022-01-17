#include "Numerical_Methods_Class.h"
#include "linspace.h"
#include "common.h"
#include "matrix_operations.h"

void Numerical_Methods_Class::RK2(int numberSpins) {

    _numberOfSpins = numberSpins;
    LinspaceClass SpinChainExchange;
    MatrixOperationsClass MatrixOperations;

    // notifies the user of what code they are running
    std::cout << "You are running the RK2 tester chainspin code.\n" << std::endl;

    std::cout << "Enter the LHS spin position for the driving region: ";
    std::cin >> _drivingRegionLHS;

    std::cout << "Enter the _stepsize: ";
    std::cin >> _stepsize;

    std::cout << "Enter the maximum number of iterations: ";
    std::cin >> _stopIterationValue; // itermax can be inputted in scientific notation or as a float
    std::cout << "\n";

    _maxSimulatedTime = _stepsize * _stopIterationValue;
    std::cout << "This will simulate a time of " << _maxSimulatedTime << "[s]." << std::endl;

    std::cout << "Enter the filename identifier: ";
    std::cin >> _fileName;
    std::cout << "\n";

    if(_drivingRegionRHS > _numberOfSpins) {
        std::cout << "The width of the domain takes it past the maximum number of spins. Exiting...";
        exit(5);
    }

    const int c_numberOfSpinPairs = _numberOfSpins - 1; // Used to indicate array lengths and tidy notation

    SpinChainExchange.set_values(_exchangeMinimum, _exchangeMaximum, c_numberOfSpinPairs,true);
    SpinChainExchange.generate_array();
    _fullChainExchangeValues = SpinChainExchange.build_spinchain();

    //TODO Turn the initial conditions lines into a separate function
    //Temporary vectors to hold the initial conditions (InitCond) of the chain along each axis. Declared separately to allow for non-isotropic conditions
    std::vector<double> mXInitCond(_numberOfSpins, _initialMagMomentX), mYInitCond(_numberOfSpins, _initialMagMomentY), mZInitCond(_numberOfSpins, _initialMagMomentZ);

    // Appends initial conditions to the vectors
    _mXEstStart.insert(_mXEstStart.end(), mXInitCond.begin(), mXInitCond.end());
    _mYEstStart.insert(_mYEstStart.end(), mYInitCond.begin(), mYInitCond.end());
    _mZEstStart.insert(_mZEstStart.end(), mZInitCond.begin(), mZInitCond.end());

    // Delete temporary vectors that held the initial conditions
    mXInitCond.clear();
    mYInitCond.clear();
    mZInitCond.clear();

    // This zero is the (N+1)th spin on the RHS of the chain
    _mXEstStart.push_back(0);
    _mYEstStart.push_back(0);
    _mZEstStart.push_back(0);

    // Creates files to save the data. All files will have (namefile) in them to make them clearly identifiable.
    std::ofstream mXFile("./Outputs_all/Output_RK2/rk2_mx_"+_fileName+".csv");
    std::ofstream mYFile("./Outputs_all/Output_RK2/rk2_my_"+_fileName+".csv");
    std::ofstream mZFile("./Outputs_all/Output_RK2/rk2_mz_"+_fileName+".csv");

    /* An increment of any RK method (such as RK4 which has k1, k2, k3 & k4) will be referred to as a stage to remove
     * confusion with the stepsize (h) which is referred to as a step or halfstep (h/2)*/
    for (long iterationIndex = _startIterationValue; iterationIndex <= (long) _stopIterationValue; iterationIndex++) {

        _totalTime += _stepsize;
        double t0 = _totalTime; // The initial time of the iteration, and the time at the first stage; the start of the interval and the first step of RK2. Often called 't0' in literature
        double t0HalfStep = _totalTime + _stepsizeHalf; // The time at the midpoint of the interval; where the second stage of RK2 is. Often called 't0*h/2' in literature

        std::vector<double> mXEstMid(_numberOfSpins+2, 0); // The estimate of the slope for the x-axis magnetic moment component at the midpoint
        std::vector<double> mYEstMid(_numberOfSpins+2, 0); // The estimate of the slope for the y-axis magnetic moment component at the midpoint
        std::vector<double> mZEstMid(_numberOfSpins+2, 0); // The estimate of the slope for the z-axis magnetic moment component at the midpoint

        for (int spin = 1; spin <= _numberOfSpins+1; spin++) { //skips the 0th and final spins as they will always be zero valued
            int LHS_spin = spin - 1, RHS_spin = spin + 1;

            // The k1 stage will use the real time (also known as t0 in https://lpsa.swarthmore.edu/NumInt/NumIntSecond.html)
            double mx_star = _mXEstStart[spin], mx_star_LHS = _mXEstStart[LHS_spin], mx_star_RHS = _mXEstStart[RHS_spin];
            double my_star = _mYEstStart[spin], my_star_LHS = _mYEstStart[LHS_spin], my_star_RHS = _mYEstStart[RHS_spin];
            double mz_star = _mZEstStart[spin], mz_star_LHS = _mZEstStart[LHS_spin], mz_star_RHS = _mZEstStart[RHS_spin];

            double k1_mx, k1_my, k1_mz;
            double Heff_x_k1, Heff_y_k1, Heff_z_k1;

            //Handles all the H_eff values. The pulse will only be on Heff_x, and it will also only be when the code is at the selected starting spin
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) { //add drive over a region by setting (spin == [lower_lim, upper_lim] for region width)
                Heff_x_k1 = _fullChainExchangeValues[LHS_spin] * mx_star_LHS + _fullChainExchangeValues[spin] * mx_star_RHS + _biasFieldDrivingAmplitude*cos(_drivingAngularFrequency * t0);
            } else {
                Heff_x_k1 = _fullChainExchangeValues[LHS_spin] * mx_star_LHS + _fullChainExchangeValues[spin] * mx_star_RHS;
            }
            Heff_y_k1 = _fullChainExchangeValues[LHS_spin] * my_star_LHS + _fullChainExchangeValues[spin] * my_star_RHS;
            Heff_z_k1 = _fullChainExchangeValues[LHS_spin] * mz_star_LHS + _fullChainExchangeValues[spin] * mz_star_RHS + _biasField;

            k1_mx = -1 * _gyroscopicMagneticConstant * (my_star * Heff_z_k1 - mz_star * Heff_y_k1);
            k1_my = +1 * _gyroscopicMagneticConstant * (mx_star * Heff_z_k1 - mz_star * Heff_x_k1);
            k1_mz = -1 * _gyroscopicMagneticConstant * (mx_star * Heff_y_k1 - my_star * Heff_x_k1);

            mXEstMid[spin] = mx_star + k1_mx*_stepsizeHalf;
            my_1[spin] = my_star + k1_my*_stepsizeHalf;
            mz_1[spin] = mz_star + k1_mz*_stepsizeHalf;
        }

        std::vector<double> mx_star_t0h(_numberOfSpins+2,0);
        std::vector<double> my_star_t0h(_numberOfSpins+2,0);
        std::vector<double> mz_star_t0h(_numberOfSpins+2,0);

        for (int spin = 1; spin <= _numberOfSpins+1; spin++) {
            int LHS_spin = spin - 1, RHS_spin = spin + 1;
            double mx1 = mx_1[spin], mx1_LHS = mx_1[LHS_spin], mx1_RHS = mx_1[RHS_spin];
            double my1 = my_1[spin], my1_LHS = my_1[LHS_spin], my1_RHS = my_1[RHS_spin];
            double mz1 = mz_1[spin], mz1_LHS = mz_1[LHS_spin], mz1_RHS = mz_1[RHS_spin];

            double k2_mx, k2_my, k2_mz;
            double Heff_x_k2, Heff_y_k2, Heff_z_k2;
            //none of the Heff values appear to be using mx,my or mz are h/2 (besides the time component)
            if (spin >= _drivingRegionLHS && spin <= _drivingRegionRHS) {
                Heff_x_k2 = _fullChainExchangeValues[LHS_spin] * mx1_LHS + _fullChainExchangeValues[spin] * mx1_RHS + _biasFieldDrivingAmplitude*cos(_drivingAngularFrequency * t0HalfStep);
            } else {
                Heff_x_k2 = _fullChainExchangeValues[LHS_spin] * mx1_LHS + _fullChainExchangeValues[spin] * mx1_RHS;
            }
            Heff_y_k2 = _fullChainExchangeValues[LHS_spin] * my1_LHS + _fullChainExchangeValues[spin] * my1_RHS;
            Heff_z_k2 = _fullChainExchangeValues[LHS_spin] * mz1_LHS + _fullChainExchangeValues[spin] * mz1_RHS + _biasField;

            k2_mx = -1 * _gyroscopicMagneticConstant * ( my1*Heff_z_k2 - mz1*Heff_y_k2 );
            k2_my = +1 * _gyroscopicMagneticConstant * ( mx1*Heff_z_k2 - mz1*Heff_x_k2 );
            k2_mz = -1 * _gyroscopicMagneticConstant * ( mx1*Heff_y_k2 - my1*Heff_x_k2 );

            //Takes above values and substitutes them into the RK4 formula: y(x+h) =
            //removed the *_stepsize from the end of each line below
            mx_star_t0h[spin] = _mXEstStart[spin] + k2_mx*_stepsize;
            my_star_t0h[spin] = _mYEstStart[spin] + k2_my*_stepsize;
            mz_star_t0h[spin] = _mZEstStart[spin] + k2_mz*_stepsize;

            //This section is to be used for debugging. It will print when a component of the magnetisation start tending to inf

            if (mx_star_t0h[spin] >= 5000){
                std::cout << "Error. Value of mx was greater than 5000 at spin(" << spin << "), iter(" << iterationIndex << ")." << std::flush;
                std::cout << " Test info as follows: _numberOfSpins = " << _numberOfSpins << "; starting spin = " << _drivingRegionLHS << "; itermax = " << _stopIterationValue << "; _stepsize: " << _stepsize << std::endl;
                exit(2);
            }

            if (my_star_t0h[spin] >= 5000){
                std::cout << "Error. Value of my was greater than 5000 at spin(" << spin << "), iter(" << iterationIndex << ")." << std::flush;
                std::cout << " Test info as follows: _numberOfSpins = " << _numberOfSpins << "; starting spin = " << _drivingRegionLHS << "; itermax = " << _stopIterationValue << "; _stepsize: " << _stepsize << std::endl;
                exit(3);
            }

            if (mz_star_t0h[spin] >= 5000){
                std::cout << "Error. Value of mz was greater than 5000 at spin(" << spin << "), iter(" << iterationIndex << ")." << std::flush;
                std::cout << " Test info as follows: _numberOfSpins = " << _numberOfSpins << "; starting spin = " << _drivingRegionLHS << "; itermax = " << _stopIterationValue << "; _stepsize: " << _stepsize << std::endl;
                exit(4);
            }

        }

        mXEstMid.clear();
        my_1.clear();
        mz_1.clear();

        _mXEstStart.clear();
        _mYEstStart.clear();
        _mZEstStart.clear();

        //all spin positions have now been calculated
        if ( iterationIndex % int(_stopIterationValue*0.01) == 0 ) {
            std::cout << "Reporting at: " << iterationIndex << std::endl;
            MatrixOperations.PrintVector(mx_star_t0h);
            for (int j = 1; j < _numberOfSpins + 1; j++) {
                mXFile << mx_star_t0h[j] << ",";
                mYFile << my_star_t0h[j] << ",";
                mZFile << mz_star_t0h[j] << ",";
                if (j == _numberOfSpins) {
                    mXFile << mx_star_t0h[j] << std::flush;
                    mYFile << my_star_t0h[j] << std::flush;
                    mZFile << mz_star_t0h[j] << std::flush;
                }
            }
            mXFile << std::endl;
            mYFile << std::endl;
            mZFile << std::endl;
        }

        _mXEstStart = mx_star_t0h;
        _mYEstStart = my_star_t0h;
        _mZEstStart = mz_star_t0h;

    }
    mXFile.close();
    mYFile.close();
    mZFile.close();

    std::cout << "Finished with: _stepsize = " << _stepsize << "; itermax = " << _stopIterationValue << "; filename = " << _fileName <<  std::endl;

}

void Numerical_Methods_Class::StreamToString() {
    // Create an output string stream
    std::ostringstream stepsizeObj;
    std::ostringstream stopiterObj;
    stepsizeObj << _stepsize;
    stopiterObj << _stopIterationValue;
    // Get string from output string stream
    _stepsizeString = stepsizeObj.str();
    std::string strObj2 = stopiterObj.str();
    stepsizeObj.clear();
    stopiterObj.clear();
}