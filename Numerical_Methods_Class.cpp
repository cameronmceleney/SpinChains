#include "Numerical_Methods_Class.h"
#include "linspace.h"
#include "common.h"
#include "matrix_operations.h"

void Numerical_Methods_Class::RK2(int numberSpins) {

    _numberOfSpins = numberSpins;
    LinspaceClass Linspace;
    MatrixOperationsClass MatrixOperations;
    
    double mx0 = 0, my0 = 0, mz0 = _magnetisationSaturation, time = _initalTime;

    // notifies the user of what code they are running
    std::cout << "You are running the RK2 tester chainspin code.\n" << std::endl;
    
    // the driving region is hardcoded in as 200 spins wide, which is why only the LHS position is requested
    std::cout << "Enter the LHS spin position for the driving region: ";
    std::cin >> _spinStart;

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

    if(_drivingRegionRHSSpin > _numberOfSpins){
        std::cout << "The width of the domain takes it past the maximum number of spins. Exiting...";
        exit(5);
    }

    const int num_spinpairs = _numberOfSpins - 1; //finds the number of pairs of spins in the chain
    Linspace.set_values(_exchangeMinimum, _exchangeMaximum, num_spinpairs,true);
    _linspaceExchangeValues = Linspace.generate_array();
    // Insert is faster for large values of numbers compared to push_back().
    std::vector<double> fullChainExchangeValues{0}; // Initialised with a zero to account for the exchange from the (P-1)th LHS spin
    fullChainExchangeValues.insert(fullChainExchangeValues.end(), _linspaceExchangeValues.begin(), _linspaceExchangeValues.end());
    fullChainExchangeValues.push_back(0); // Appends a zero to the end to account for the exchange from the (N+1)th RHS spin

    // Creates two sets of vectors. The _prev vectors are needed for the program and represent the values of the
    // previous iteration of calculations. The init vectors are temporary vectors used to contain the values declared
    // in the initial state of the program. The zero that the _prev vectors contain is the (P-1)th spin
    std::vector<double> mx_star_t0{0}, my_star_t0{0}, mz_star_t0{0}, init_mx(_numberOfSpins, mx0), init_my(_numberOfSpins, my0), init_mz(_numberOfSpins, mz0);
    //appends initial conditions to the vectors
    mx_star_t0.insert(mx_star_t0.end(), init_mx.begin(), init_mx.end());
    my_star_t0.insert(my_star_t0.end(), init_my.begin(), init_my.end());
    mz_star_t0.insert(mz_star_t0.end(), init_mz.begin(), init_mz.end());

    //Add zeros to the end of the vectors. This zero is the N+1th spin
    mx_star_t0.push_back(0);
    my_star_t0.push_back(0);
    mz_star_t0.push_back(0);

    //clears the temporary vectors that housed the initial values
    init_mx.clear();
    init_my.clear();
    init_mz.clear();

    // Create an output string stream
    std::ostringstream streamObj;
    std::ostringstream streamObj2;
    //Add doubles to stream
    streamObj << _stepsize;
    streamObj2 << _stopIterationValue;
    // Get string from output string stream
    std::string strObj = streamObj.str();
    std::string strObj2 = streamObj2.str();
    streamObj.clear();
    streamObj2.clear();

    // std::string filename_iden = "_iter-"+strObj2+"_step-"+strObj;

    // Creates files to save the data. All files will have (namefile) in them to make them clearly identifiable.
    std::ofstream mxfile("./Outputs_all/Output_RK2/rk2_mx_"+_fileName+".csv");
    std::ofstream myfile("./Outputs_all/Output_RK2/rk2_my_"+_fileName+".csv");
    std::ofstream mzfile("./Outputs_all/Output_RK2/rk2_mz_"+_fileName+".csv");

    //debugging section before loop starts
    //for (int s=0; s < num_spinpairs+2; s++) {std::cout << "Elem(" << s << ")=" << fullChainExchangeValues[s] << ",\t" << std::flush;} std::cout <<'\n';
    //for (int t=0; t < _numberOfSpins+2; t++) {std::cout << "Elem(" << t << ")=" << mz_prev[t] << ",\t";} std::cout <<'\n';
    //std::cout << "Size Mx: " << mx_prev.size() << "| Size JI_Values: " << fullChainExchangeValues.size() << std::endl;

    for (long iteration = _startIterationValue; iteration <= (long) _stopIterationValue; iteration++) {

        //Keeps track of the real time. This variable is used for anything that varying in time, like the driving pulses
        time += _stepsize;
        double t0 = time, half_step = _stepsize/2, t0half = time + half_step;

        std::vector<double> mx_1(_numberOfSpins+2, 0);
        std::vector<double> my_1(_numberOfSpins+2, 0);
        std::vector<double> mz_1(_numberOfSpins+2, 0);

        for (int spin = 1; spin <= _numberOfSpins+1; spin++) { //skips the 0th and final spins as they will always be zero valued
            int LHS_spin = spin - 1, RHS_spin = spin + 1;

            // The k1 stage will use the real time (also known as t0 in https://lpsa.swarthmore.edu/NumInt/NumIntSecond.html##section6)
            double mx_star = mx_star_t0[spin], mx_star_LHS = mx_star_t0[LHS_spin], mx_star_RHS = mx_star_t0[RHS_spin];
            double my_star = my_star_t0[spin], my_star_LHS = my_star_t0[LHS_spin], my_star_RHS = my_star_t0[RHS_spin];
            double mz_star = mz_star_t0[spin], mz_star_LHS = mz_star_t0[LHS_spin], mz_star_RHS = mz_star_t0[RHS_spin];

            double k1_mx, k1_my, k1_mz;
            double Heff_x_k1, Heff_y_k1, Heff_z_k1;

            //Handles all the H_eff values. The pulse will only be on Heff_x, and it will also only be when the code is at the selected starting spin
            if (spin >= _drivingRegionLHSSpin && spin <= _drivingRegionRHSSpin) { //add drive over a region by setting (spin == [lower_lim, upper_lim] for region width)
                Heff_x_k1 = fullChainExchangeValues[LHS_spin] * mx_star_LHS + fullChainExchangeValues[spin] * mx_star_RHS + _biasFieldDrivingAmplitude*cos(_drivingAngularFrequency * t0);
            } else {
                Heff_x_k1 = fullChainExchangeValues[LHS_spin] * mx_star_LHS + fullChainExchangeValues[spin] * mx_star_RHS;
            }
            Heff_y_k1 = fullChainExchangeValues[LHS_spin] * my_star_LHS + fullChainExchangeValues[spin] * my_star_RHS;
            Heff_z_k1 = fullChainExchangeValues[LHS_spin] * mz_star_LHS + fullChainExchangeValues[spin] * mz_star_RHS + _biasField;

            k1_mx = -1 * _gyroscopicMagneticConstant * (my_star * Heff_z_k1 - mz_star * Heff_y_k1);
            k1_my = +1 * _gyroscopicMagneticConstant * (mx_star * Heff_z_k1 - mz_star * Heff_x_k1);
            k1_mz = -1 * _gyroscopicMagneticConstant * (mx_star * Heff_y_k1 - my_star * Heff_x_k1);

            mx_1[spin] = mx_star + k1_mx*half_step;
            my_1[spin] = my_star + k1_my*half_step;
            mz_1[spin] = mz_star + k1_mz*half_step;
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
            if (spin >= _drivingRegionLHSSpin && spin <= _drivingRegionRHSSpin) {
                Heff_x_k2 = fullChainExchangeValues[LHS_spin] * mx1_LHS + fullChainExchangeValues[spin] * mx1_RHS + _biasFieldDrivingAmplitude*cos(_drivingAngularFrequency * t0half);
            } else {
                Heff_x_k2 = fullChainExchangeValues[LHS_spin] * mx1_LHS + fullChainExchangeValues[spin] * mx1_RHS;
            }
            Heff_y_k2 = fullChainExchangeValues[LHS_spin] * my1_LHS + fullChainExchangeValues[spin] * my1_RHS;
            Heff_z_k2 = fullChainExchangeValues[LHS_spin] * mz1_LHS + fullChainExchangeValues[spin] * mz1_RHS + _biasField;

            k2_mx = -1 * _gyroscopicMagneticConstant * ( my1*Heff_z_k2 - mz1*Heff_y_k2 );
            k2_my = +1 * _gyroscopicMagneticConstant * ( mx1*Heff_z_k2 - mz1*Heff_x_k2 );
            k2_mz = -1 * _gyroscopicMagneticConstant * ( mx1*Heff_y_k2 - my1*Heff_x_k2 );

            //Takes above values and substitutes them into the RK4 formula: y(x+h) =
            //removed the *_stepsize from the end of each line below
            mx_star_t0h[spin] = mx_star_t0[spin] + k2_mx*_stepsize;
            my_star_t0h[spin] = my_star_t0[spin] + k2_my*_stepsize;
            mz_star_t0h[spin] = mz_star_t0[spin] + k2_mz*_stepsize;

            //This section is to be used for debugging. It will print when a component of the magnetisation start tending to inf

            if (mx_star_t0h[spin] >= 5000){
                std::cout << "Error. Value of mx was greater than 5000 at spin(" << spin << "), iter(" << iteration << ")." << std::flush;
                std::cout << " Test info as follows: _numberOfSpins = " << _numberOfSpins << "; starting spin = " << _spinStart << "; itermax = " << _stopIterationValue << "; _stepsize: " << _stepsize << std::endl;
                exit(2);
            }

            if (my_star_t0h[spin] >= 5000){
                std::cout << "Error. Value of my was greater than 5000 at spin(" << spin << "), iter(" << iteration << ")." << std::flush;
                std::cout << " Test info as follows: _numberOfSpins = " << _numberOfSpins << "; starting spin = " << _spinStart << "; itermax = " << _stopIterationValue << "; _stepsize: " << _stepsize << std::endl;
                exit(3);
            }

            if (mz_star_t0h[spin] >= 5000){
                std::cout << "Error. Value of mz was greater than 5000 at spin(" << spin << "), iter(" << iteration << ")." << std::flush;
                std::cout << " Test info as follows: _numberOfSpins = " << _numberOfSpins << "; starting spin = " << _spinStart << "; itermax = " << _stopIterationValue << "; _stepsize: " << _stepsize << std::endl;
                exit(4);
            }

        }

        mx_1.clear();
        my_1.clear();
        mz_1.clear();

        mx_star_t0.clear();
        my_star_t0.clear();
        mz_star_t0.clear();

        //all spin positions have now been calculated
        if ( iteration % int(_stopIterationValue*0.01) == 0 ) {
            std::cout << "Reporting at: " << iteration << std::endl;
            MatrixOperations.PrintVector(mx_star_t0h);
            for (int j = 1; j < _numberOfSpins + 1; j++) {
                mxfile << mx_star_t0h[j] << ",";
                myfile << my_star_t0h[j] << ",";
                mzfile << mz_star_t0h[j] << ",";
                if (j == _numberOfSpins) {
                    mxfile << mx_star_t0h[j] << std::flush;
                    myfile << my_star_t0h[j] << std::flush;
                    mzfile << mz_star_t0h[j] << std::flush;
                }
            }
            mxfile << std::endl;
            myfile << std::endl;
            mzfile << std::endl;
        }

        mx_star_t0 = mx_star_t0h;
        my_star_t0 = my_star_t0h;
        mz_star_t0 = mz_star_t0h;

    }
    mxfile.close();
    myfile.close();
    mzfile.close();

    std::cout << "Finished with: _stepsize = " << _stepsize << "; itermax = " << _stopIterationValue << "; filename = " << _fileName <<  std::endl;

}