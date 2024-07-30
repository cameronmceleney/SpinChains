//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/DemagnetisationFields.h"

DemagnetisationFields::DemagnetisationFields(SimulationParameters* sharedSimParams,
                                             SimulationStates* sharedSimStates,
                                             SimulationFlags* sharedSimFlags)

   : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}

void DemagnetisationFields::initialise() {

    // Need length, width, and thickness. For now will just create these as literals (in m to allow safe use of float)
    auto length = static_cast<float>(_simParams->systemTotalSpins * _simParams->latticeConstant * 1e9);
    float width = 2.0f;
    float thickness = 8.0f;

    if (_simStates->demagFactors.isDemagFactorDefaultState()) {
        _simStates->demagFactors.Nz = _calculateDemagFactorsUniformPrism(length, width, thickness);
        _simStates->demagFactors.Ny = _calculateDemagFactorsUniformPrism(thickness, length, width);
        _simStates->demagFactors.Nx = _calculateDemagFactorsUniformPrism(width, thickness, length);
    }
    // Check if the demagnetisation factors are acceptable
    if (!_simStates->demagFactors.isDemagFactorErrorAcceptable()) {
        std::cerr << "Error: Demagnetisation factors are not acceptable." << std::endl;
        std::exit(1); // or consider throwing an exception
    }

    // std::cout << _simStates->demagFactors.Nx << ", " << _simStates->demagFactors.Ny << ", " << _simStates->demagFactors.Nz << std::endl;
}

void DemagnetisationFields::DemagnetisationFieldIntense( const std::vector<double> &mxTerms,
                                                         const std::vector<double> &myTerms,
                                                         const std::vector<double> &mzTerms, std::vector<double> &H_dx,
                                                         std::vector<double> &H_dy, std::vector<double> &H_dz ) {

    std::vector<double> totalHd_x(GV.GetNumSpins() + 2, 0), totalHd_y(GV.GetNumSpins() + 2, 0), totalHd_z(GV.GetNumSpins() + 2, 0);

    // Loop over all sites
    for(int currentSite = 1; currentSite <= GV.GetNumSpins(); currentSite++) {
        double localHd_x = 0.0;
        double localHd_y = 0.0;
        double localHd_z = 0.0;

        // Convolution-like sum over all other sites
        for(int otherSites = 1; otherSites <= GV.GetNumSpins(); otherSites++) {
            if(currentSite != otherSites) { // Avoid self-interaction
                localHd_x += _simStates->demagFactors.Nx * (mxTerms[currentSite] - mxTerms[otherSites]);
                localHd_y += _simStates->demagFactors.Ny * (myTerms[currentSite] - myTerms[otherSites]);
                localHd_z += _simStates->demagFactors.Nz * (mzTerms[currentSite] - mzTerms[otherSites]);
            }
        }

        totalHd_x[currentSite] += localHd_x;
        totalHd_y[currentSite] += localHd_y;
        totalHd_z[currentSite] += localHd_z;
    }

    /*
     * Do not average over the number of spins! Here the H_d components for every site are stored as the total H_d
     * experienced by each site. Averaging is therefore unnecessary unless only a single totalH_d is returned by the
     * function and it is the sum of every demag component in the system (1/n SUM_i H_d_x_i).
     */

    // Keeping separate for now to aid debugging!
    H_dx = totalHd_x;
    H_dy = totalHd_y;
    H_dz = totalHd_z;
}

void DemagnetisationFields::DemagField1DComplex(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                           std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                           int iteration, std::string rkStageName) {

    int gotNumSpins = GV.GetNumSpins();  // Vectors and arrays defined in function don't include pinned end terms; 4000 sites in length
    int trueNumSpins = GV.GetNumSpins() + 2;  // Vectors and arrays defined out with function include pinned end terms; 4002 sites in length

    const double imagTerm = 0.0, normalizationFactor = 1.0 / static_cast<double>(gotNumSpins);

    auto fftw_alloc_and_check = [](const char* var_name, const int& size) -> fftw_complex* {
        /*
         * Temporary lambda function for debugging. Keep within `DemagField1D` until debugging complete.
         * Allocates memory for FFTW-suitable arrays.Currently being cautious, so throw exception if allocation fails,
         * and else set memory to zero to ensure data integrity.
         */
        auto *mTerm = static_cast<fftw_complex*>(fftw_alloc_complex(size));
        if (!mTerm) {
            throw std::runtime_error(std::string("Failed to allocate memory for ") + var_name);
        } else {
            std::memset(mTerm, 0, sizeof(fftw_complex) * size);
            return mTerm;
        }
    };

    auto PrintFFTWVector = [](int numSpins, fftw_complex* vecToPrint, const char* vecName, double testValReal, double testValImag) {
        /*
         * Temporary lambda function for debugging. Keep within `DemagField1D` until debugging complete.
         * Prints the contents of an FFTW-suitable array. Helps find abnormal and unexpected values (compared to my
         * analytical solutions.
         */
        for (int i = 0; i < numSpins; i++) {
            if (vecToPrint[i][0] != testValReal)
                std::cout << vecName << "[" << i << "][0] = " << vecToPrint[i][0] << std::endl;
            else if (vecToPrint[i][1] != testValImag)
                std::cout << vecName << "[" << i << "][1] = " << vecToPrint[i][1] << std::endl;
        }

    };

    // Assign memory for FFTW-suitable arrays for magnetic moments; used during computation and for RMSE calculation
    auto *mX = fftw_alloc_and_check("mx", gotNumSpins);
    auto *mY = fftw_alloc_and_check("my", gotNumSpins);
    auto *mZ = fftw_alloc_and_check("mz", gotNumSpins);

    // Additional memory for demagnetisation field components; used for output; initialise at zero (empty)
    auto *hdX = fftw_alloc_and_check("hdX", gotNumSpins);
    auto *hdY = fftw_alloc_and_check("hdY", gotNumSpins);
    auto *hdZ = fftw_alloc_and_check("hdZ", gotNumSpins);

    // Population of memory for FFTW-suitable arrays
    for(int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        int scaled_len = currentSite + 1;  // mX (etc) is 4000 sites whereas inMxTerms (etc) is 4002 sites
        mX[currentSite][0] = inMxTerms[scaled_len]; mX[currentSite][1] = imagTerm; // CHECKED: all mX elements are [0, 0] as expected
        mY[currentSite][0] = inMyTerms[scaled_len]; mY[currentSite][1] = imagTerm; // CHECKED: all mY elements are [0, 0] as expected
        mZ[currentSite][0] = inMzTerms[scaled_len]; mZ[currentSite][1] = imagTerm; // CHECKED: all mZ elements are [1, 0] as expected
        // hdX (etc) should already be zeroes. Keeping this for debugging purposes to ensure data integrity; will remove when working
        // hdX[currentSite][0] = imagTerm; hdX[currentSite][1] = imagTerm; // CHECKED: all hdX elements are [0, 0] as expected
        // hdY[currentSite][0] = imagTerm; hdY[currentSite][1] = imagTerm; // CHECKED: all hdY elements are [0, 0] as expected
        // hdZ[currentSite][0] = imagTerm; hdZ[currentSite][1] = imagTerm; // CHECKED: all hdZ elements are [0, 0] as expected
    }
    // Have debugging mX (etc) and hdX (etc) and found that they are all being populated correctly


    // Create FFTW plans **after** population of memory
    fftw_plan planDFTForwardMx = fftw_plan_dft_1d(gotNumSpins, mX, mX, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan planDFTForwardMy = fftw_plan_dft_1d(gotNumSpins, mY, mY, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan planDFTForwardMz = fftw_plan_dft_1d(gotNumSpins, mZ, mZ, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan planDFTForwardHdX = fftw_plan_dft_1d(gotNumSpins, hdX, hdX, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan planDFTForwardHdY = fftw_plan_dft_1d(gotNumSpins, hdY, hdY, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan planDFTForwardHdZ = fftw_plan_dft_1d(gotNumSpins, hdZ, hdZ, FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute the plans **after** their definitions
    fftw_execute(planDFTForwardMx);
    fftw_execute(planDFTForwardMy);
    fftw_execute(planDFTForwardMz);
    fftw_execute(planDFTForwardHdX);
    fftw_execute(planDFTForwardHdY);
    fftw_execute(planDFTForwardHdZ);

    // Population in Fourier space (H_d = -N * M). Keep mX/mY/mZ un-mutated for RMSE calculation later
    for (int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        for (int part = 0; part < 2; part++) {
            hdX[currentSite][part] =  -1 * _simStates->demagFactors.Nx * mX[currentSite][part];  // CHECKED: all hdX elements are [0, 0] as expected
            hdY[currentSite][part] =  -1 * _simStates->demagFactors.Ny * mY[currentSite][part];  // CHECKED: all hdY elements are [0, 0] as expected
            hdZ[currentSite][part] =  -1 * _simStates->demagFactors.Nz * mZ[currentSite][part];  // CHECKED: all hdZ elements are [0, 0] TODO apart from hdZ[0][0] = [-2000, 0]. Why?!
        }
    }

    // Create inverse FFTW plans **after** population of memory
    fftw_plan planDFTBackMx = fftw_plan_dft_1d(gotNumSpins, mX, mX, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan planDFTBackMy = fftw_plan_dft_1d(gotNumSpins, mY, mY, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan planDFTBackMz = fftw_plan_dft_1d(gotNumSpins, mZ, mZ, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdX = fftw_plan_dft_1d(gotNumSpins, hdX, hdX, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdY = fftw_plan_dft_1d(gotNumSpins, hdY, hdY, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdZ = fftw_plan_dft_1d(gotNumSpins, hdZ, hdZ, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Execute the plans **after** their definitions
    fftw_execute(planDFTBackMx);
    fftw_execute(planDFTBackMy);
    fftw_execute(planDFTBackMz);
    fftw_execute(planDFTBackHdX);
    fftw_execute(planDFTBackHdY);
    fftw_execute(planDFTBackHdZ);

    // Normalise, after FFT, all arrays that were taken into Fourier space
    for (int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        for (int part = 0; part < 1; part++) {
            mX[currentSite][part] *= normalizationFactor;  // CHECKED: all mX elements are [0, 0] as expected
            mY[currentSite][part] *= normalizationFactor; // CHECKED: all mY elements are [0, 0] as expected
            mZ[currentSite][part] *= normalizationFactor; // CHECKED: all mZ elements are [1, 0] as expected
            hdX[currentSite][part] *= normalizationFactor; // CHECKED: all hdX elements are [0, 0] as expected
            hdY[currentSite][part] *= normalizationFactor; // CHECKED: all hdY elements are [0, 0] as expected
            hdZ[currentSite][part] *= normalizationFactor; // CHECKED: all hdZ elements are [-0.5, 0] as expected
        }
    }

    // Compute Root Mean Square Error (RMSE) for each component of the magnetic moment
    double rmse_mx = 0, rmse_my = 0, rmse_mz = 0;
    for (int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        int scaled_len = currentSite + 1;  // mx.size() is 4000 whereas mxTerms.size() is  4002
        // Finding SUM_{i=1}^{N}(data_in[i] - data_reconstructed[i])**2
        rmse_mx += pow(inMxTerms[scaled_len] - mX[currentSite][0], 2);  // CHECKED: rmse_mx is equal to zero
        rmse_my += pow(inMyTerms[scaled_len] - mY[currentSite][0], 2);  // CHECKED: rmse_my is equal to zero
        rmse_mz += pow(inMzTerms[scaled_len] - mZ[currentSite][0], 2);  // CHECKED: rmse_mz is equal to zero
    }

    // Finishing RMSE calculation: RMSE = Sqrt(SUM_{i=1}^{N} * 1/N)
    rmse_mx = sqrt(rmse_mx / gotNumSpins);  // CHECKED: rmse_mx is equal to zero
    rmse_my = sqrt(rmse_my / gotNumSpins);  // CHECKED: rmse_my is equal to zero
    rmse_mz = sqrt(rmse_mz / gotNumSpins);  // CHECKED: rmse_mz is equal to zero


    // FFT error estimation based on machine epsilon and FFTW's scale factor
    const double machineEpsilon = std::numeric_limits<double>::epsilon();
    // const double fftError = machineEpsilon * gotNumSpins * sqrt(gotNumSpins);
    // std::cout << fftError << std::endl;std::exit(0);

    bool mxRMSETest = false, myRMSETest = false, mzRMSETest = false;
    if (fabs(rmse_mx) < machineEpsilon)
        mxRMSETest = true;
    if (fabs(rmse_my) < machineEpsilon)
        myRMSETest = true;
    if (fabs(rmse_mz) < machineEpsilon)
        mzRMSETest = true;

    for (int trueSite = 0; trueSite < gotNumSpins + 2; trueSite++) {
        int scale_site = trueSite - 1;
        if (trueSite == 0 || trueSite == (trueNumSpins - 1) ) {
            // Boundary conditions of output arrays must always be zero; do this to ensure data integrity
            outDemagX[trueSite] = 0.0;
            outDemagY[trueSite] = 0.0;
            outDemagZ[trueSite] = 0.0;
            continue;
        }
        /*
         * outDemagX/outDemagY/outDemagZ are 4002 in length while mX/mY/mZ are 4000 in length. mX (etc) are unmutated throughout
         * the FFT, which inMxTerms (etc) are the original values. If the absolute difference between inMxTerms and mX (etc)
         * is greater than their rmse_mx (etc), then the overwriting of outDemagX (etc) by hdX (etc) should be stopped. This
         * is because the values are likely to just be noise. This is a very crude way of doing this, but it works for now.
         *
         * Caution! outDemagX (etc) are references so overwriting their elements changes what the main program sees.
         */
        outDemagX[trueSite] = mxRMSETest ? hdX[scale_site][0] : 0.0; // hdX[scale_site][0]; //(fabs((inMxTerms[trueSite] - mX[scale_site][0])) > fabs(rmse_mx)) ? hdX[scale_site][0] : 0.0;  // HELP HERE
        outDemagY[trueSite] = myRMSETest ? hdY[scale_site][0] : 0.0;// hd  Y[scale_site][0]; //(fabs((inMyTerms[trueSite] - mY[scale_site][0])) > fabs(rmse_my)) ? hdY[scale_site][0] : 0.0;  // HELP HERE
        outDemagZ[trueSite] = mzRMSETest ? hdZ[scale_site][0] : 0.0;// hdZ[scale_site][0]; //(fabs((inMzTerms[trueSite] - mZ[scale_site][0])) > fabs(rmse_mz)) ? hdZ[scale_site][0] : 0.0;  // HELP HERE
    }
    bool testOutRMS = false;
    if (testOutRMS) {
        if (_simParams->iterationEnd >= 100 && iteration % (_simParams->iterationEnd / 1000) == 0) {
            if (mxRMSETest || myRMSETest || mzRMSETest) {
                std::cout << "Iter. #" << iteration << " | RK" << rkStageName << " ";
                if (mxRMSETest)
                    std::cout << "| RMSE. mx: " << rmse_mx << " ";
                if (myRMSETest)
                    std::cout << "| RMSE. my: " << rmse_my << " ";
                if (mzRMSETest)
                    std::cout << "| RMSE. mz: " << rmse_mz << " ";
                std::cout << std::endl;
            }
            //std::cout << "Iteration #" << iteration <<" | RMSE. mx: " << rmse_mx << " | my: " << rmse_my << " | mz:  " << rmse_mz << std::endl;  // Keep for debugging
        }
    }
    /*
    std::cout << "HERE IN DEMAGFIELD1D: X" << std::endl;
    PrintVector(outDemagX, false);

    std::cout << "HERE IN DEMAGFIELD1D: Y" << std::endl;
    PrintVector(outDemagY, false);

    std::cout << "HERE IN DEMAGFIELD1D: Z" << std::endl;
    PrintVector(outDemagZ, false);
     */

    // EASY FIND
    // Clean-up. Probably could free memory for planDFTForwardMx (etc) earlier in function, but it's safer to be here
    fftw_destroy_plan(planDFTForwardMx);
    fftw_destroy_plan(planDFTForwardMy);
    fftw_destroy_plan(planDFTForwardMz);
    fftw_destroy_plan(planDFTForwardHdX);
    fftw_destroy_plan(planDFTForwardHdY);
    fftw_destroy_plan(planDFTForwardHdZ);
    fftw_destroy_plan(planDFTBackMx);
    fftw_destroy_plan(planDFTBackMy);
    fftw_destroy_plan(planDFTBackMz);
    fftw_destroy_plan(planDFTBackHdX);
    fftw_destroy_plan(planDFTBackHdY);
    fftw_destroy_plan(planDFTBackHdZ);
    fftw_free(mX);
    fftw_free(mY);
    fftw_free(mZ);
    fftw_free(hdX);
    fftw_free(hdY);
    fftw_free(hdZ);
}

void DemagnetisationFields::DemagField1DReal(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                           std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                           int iteration, std::string rkStageName) {

    int gotNumSpins = GV.GetNumSpins();  // Vectors and arrays defined in function don't include pinned end terms; 4000 sites in length
    int trueNumSpins = GV.GetNumSpins() + 2;  // Vectors and arrays defined out with function include pinned end terms; 4002 sites in length
    const double imagTerm = 0.0, normalizationFactor = 1.0 / static_cast<double>(gotNumSpins);

    auto fftw_alloc_and_check = [](const char* var_name, const int& size) -> double* {
        /*
         * Temporary lambda function for debugging. Keep within `DemagField1D` until debugging complete.
         * Allocates memory for FFTW-suitable arrays.Currently being cautious, so throw exception if allocation fails,
         * and else set memory to zero to ensure data integrity.
         */
        auto *mTerm = static_cast<double*>(fftw_alloc_real(size));
        if (!mTerm) {
            throw std::runtime_error(std::string("Failed to allocate memory for ") + var_name);
        } else {
            std::memset(mTerm, 0, sizeof(double) * size);
            return mTerm;
        }
    };

    auto PrintFFTWVector = [](int numSpins, fftw_complex* vecToPrint, const char* vecName, double testValReal, double testValImag) {
        /*
         * Temporary lambda function for debugging. Keep within `DemagField1D` until debugging complete.
         * Prints the contents of an FFTW-suitable array. Helps find abnormal and unexpected values (compared to my
         * analytical solutions.
         */
        for (int i = 0; i < numSpins; i++) {
            if (vecToPrint[i][0] != testValReal)
                std::cout << vecName << "[" << i << "][0] = " << vecToPrint[i][0] << std::endl;
            else if (vecToPrint[i][1] != testValImag)
                std::cout << vecName << "[" << i << "][1] = " << vecToPrint[i][1] << std::endl;
        }

    };

    // Assign memory for FFTW-suitable arrays for magnetic moments; used during computation and for RMSE calculation
    auto *mX = fftw_alloc_and_check("mx", gotNumSpins);
    auto *mY = fftw_alloc_and_check("my", gotNumSpins);
    auto *mZ = fftw_alloc_and_check("mz", gotNumSpins);

    // Additional memory for demagnetisation field components; used for output; initialise at zero (empty)
    auto *hdX = fftw_alloc_and_check("hdX", gotNumSpins);
    auto *hdY = fftw_alloc_and_check("hdY", gotNumSpins);
    auto *hdZ = fftw_alloc_and_check("hdZ", gotNumSpins);

    // Population of memory for FFTW-suitable arrays
    for(int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        int scaled_len = currentSite + 1;  // mX (etc) is 4000 sites whereas inMxTerms (etc) is 4002 sites
        mX[currentSite] = inMxTerms[scaled_len]; // CHECKED: all mX elements are [0, 0] as expected
        mY[currentSite] = inMyTerms[scaled_len]; // CHECKED: all mY elements are [0, 0] as expected
        mZ[currentSite] = inMzTerms[scaled_len]; // CHECKED: all mZ elements are [1, 0] as expected
        // hdX (etc) should already be zeroes. Kee; to ensure data integrity; will remove when working
        // hdX[currentSite][0] = imagTerm; hdX[currentSite][1] = imagTerm; // CHECKED: all hdX elements are [0, 0] as expected
        // hdY[currentSite][0] = imagTerm; hdY[currentSite][1] = imagTerm; // CHECKED: all hdY elements are [0, 0] as expected
        // hdZ[currentSite][0] = imagTerm; hdZ[currentSite][1] = imagTerm; // CHECKED: all hdZ elements are [0, 0] as expected
    }
    // Have debugging mX (etc) and hdX (etc) and found that they are all being populated correctly


    // Create FFTW plans **after** population of memory
    fftw_plan planDFTForwardMx = fftw_plan_r2r_1d(gotNumSpins, mX, mX, FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan planDFTForwardMy = fftw_plan_r2r_1d(gotNumSpins, mY, mY, FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan planDFTForwardMz = fftw_plan_r2r_1d(gotNumSpins, mZ, mZ, FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan planDFTForwardHdX = fftw_plan_r2r_1d(gotNumSpins, hdX, hdX, FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan planDFTForwardHdY = fftw_plan_r2r_1d(gotNumSpins, hdY, hdY, FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan planDFTForwardHdZ = fftw_plan_r2r_1d(gotNumSpins, hdZ, hdZ, FFTW_R2HC, FFTW_ESTIMATE);

    // Execute the plans **after** their definitions
    fftw_execute(planDFTForwardMx);
    fftw_execute(planDFTForwardMy);
    fftw_execute(planDFTForwardMz);
    fftw_execute(planDFTForwardHdX);
    fftw_execute(planDFTForwardHdY);
    fftw_execute(planDFTForwardHdZ);

    // Population in Fourier space (H_d = -N * M). Keep mX/mY/mZ un-mutated for RMSE calculation later
    for (int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        hdX[currentSite] = -1 * _simStates->demagFactors.Nx * mX[currentSite];  // CHECKED: all hdX elements are [0, 0] as expected
        hdY[currentSite] = -1 * _simStates->demagFactors.Ny * mY[currentSite];  // CHECKED: all hdY elements are [0, 0] as expected
        hdZ[currentSite] = -1 * _simStates->demagFactors.Nz * mZ[currentSite];  // CHECKED: all hdZ elements are [0, 0] TODO apart from hdZ[0][0] = [-2000, 0]. Why?!

    }

    // Create inverse FFTW plans **after** population of memory
    fftw_plan planDFTBackMx = fftw_plan_r2r_1d(gotNumSpins, mX, mX, FFTW_HC2R, FFTW_ESTIMATE);
    fftw_plan planDFTBackMy = fftw_plan_r2r_1d(gotNumSpins, mY, mY, FFTW_HC2R, FFTW_ESTIMATE);
    fftw_plan planDFTBackMz = fftw_plan_r2r_1d(gotNumSpins, mZ, mZ, FFTW_HC2R, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdX = fftw_plan_r2r_1d(gotNumSpins, hdX, hdX, FFTW_HC2R, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdY = fftw_plan_r2r_1d(gotNumSpins, hdY, hdY, FFTW_HC2R, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdZ = fftw_plan_r2r_1d(gotNumSpins, hdZ, hdZ, FFTW_HC2R, FFTW_ESTIMATE);

    // Execute the plans **after** their definitions
    fftw_execute(planDFTBackMx);
    fftw_execute(planDFTBackMy);
    fftw_execute(planDFTBackMz);
    fftw_execute(planDFTBackHdX);
    fftw_execute(planDFTBackHdY);
    fftw_execute(planDFTBackHdZ);

    // Normalise, after FFT, all arrays that were taken into Fourier space
    for (int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        mX[currentSite] *= normalizationFactor;  // CHECKED: all mX elements are [0, 0] as expected
        mY[currentSite] *= normalizationFactor; // CHECKED: all mY elements are [0, 0] as expected
        mZ[currentSite] *= normalizationFactor; // CHECKED: all mZ elements are [1, 0] as expected
        hdX[currentSite] *= normalizationFactor; // CHECKED: all hdX elements are [0, 0] as expected
        hdY[currentSite] *= normalizationFactor; // CHECKED: all hdY elements are [0, 0] as expected
        hdZ[currentSite] *= normalizationFactor; // CHECKED: all hdZ elements are [-0.5, 0] as expected
    }

    // Compute Root Mean Square Error (RMSE) for each component of the magnetic moment
    double rmse_mx = 0, rmse_my = 0, rmse_mz = 0;
    for (int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        int scaled_len = currentSite + 1;  // mx.size() is 4000 whereas mxTerms.size() is  4002
        // Finding SUM_{i=1}^{N}(data_in[i] - data_reconstructed[i])**2
        rmse_mx += pow(inMxTerms[scaled_len] - mX[currentSite], 2);  // CHECKED: rmse_mx is equal to zero
        rmse_my += pow(inMyTerms[scaled_len] - mY[currentSite], 2);  // CHECKED: rmse_my is equal to zero
        rmse_mz += pow(inMzTerms[scaled_len] - mZ[currentSite], 2);  // CHECKED: rmse_mz is equal to zero
    }

    // Finishing RMSE calculation: RMSE = Sqrt(SUM_{i=1}^{N} * 1/N)
    rmse_mx = sqrt(rmse_mx / gotNumSpins);  // CHECKED: rmse_mx is equal to zero
    rmse_my = sqrt(rmse_my / gotNumSpins);  // CHECKED: rmse_my is equal to zero
    rmse_mz = sqrt(rmse_mz / gotNumSpins);  // CHECKED: rmse_mz is equal to zero


    // FFT error estimation based on machine epsilon and FFTW's scale factor
    const double machineEpsilon = std::numeric_limits<double>::epsilon();
    // const double fftError = machineEpsilon * gotNumSpins * sqrt(gotNumSpins);
    // std::cout << fftError << std::endl;std::exit(0);

    bool mxRMSETest, myRMSETest, mzRMSETest;
    if (fabs(rmse_mx) < machineEpsilon) { mxRMSETest = true; }
    else {mxRMSETest = false;}
    if (fabs(rmse_my) < machineEpsilon) { myRMSETest = true; }
    else { myRMSETest = false; }
    if (fabs(rmse_mz) < machineEpsilon) { mzRMSETest = true; }
    else { mzRMSETest = false; }

    double combinedRMSE = sqrt((pow(rmse_mx, 2) + pow(rmse_my, 2) + pow(rmse_mz, 2)) / 3.0 );
    bool applyDemag;
    if (combinedRMSE < machineEpsilon) {applyDemag = true;}
    else {applyDemag = false;}
    // if (mxRMSETest || myRMSETest || mzRMSETest) {applyDemag = true;}
    // else {applyDemag = false;}

    for (int trueSite = 0; trueSite < gotNumSpins + 2; trueSite++) {
        int scale_site = trueSite - 1;
        if (trueSite == 0 || trueSite == (trueNumSpins - 1) ) {
            // Boundary conditions of output arrays must always be zero; do this to ensure data integrity
            outDemagX[trueSite] = 0.0;
            outDemagY[trueSite] = 0.0;
            outDemagZ[trueSite] = 0.0;
            continue;
        }
        /*
         * outDemagX/outDemagY/outDemagZ are 4002 in length while mX/mY/mZ are 4000 in length. mX (etc) are unmutated throughout
         * the FFT, which inMxTerms (etc) are the original values. If the absolute difference between inMxTerms and mX (etc)
         * is greater than their rmse_mx (etc), then the overwriting of outDemagX (etc) by hdX (etc) should be stopped. This
         * is because the values are likely to just be noise. This is a very crude way of doing this, but it works for now.
         *
         * Caution! outDemagX (etc) are references so overwriting their elements changes what the main program sees.
         */
        outDemagX[trueSite] = applyDemag ? hdX[scale_site] : 0.0; // hdX[scale_site][0]; //(fabs((inMxTerms[trueSite] - mX[scale_site][0])) > fabs(rmse_mx)) ? hdX[scale_site][0] : 0.0;  // HELP HERE
        outDemagY[trueSite] = applyDemag ? hdY[scale_site] : 0.0;// hdY[scale_site][0]; //(fabs((inMyTerms[trueSite] - mY[scale_site][0])) > fabs(rmse_my)) ? hdY[scale_site][0] : 0.0;  // HELP HERE
        outDemagZ[trueSite] = applyDemag ? hdZ[scale_site] : 0.0;// hdZ[scale_site][0]; //(fabs((inMzTerms[trueSite] - mZ[scale_site][0])) > fabs(rmse_mz)) ? hdZ[scale_site][0] : 0.0;  // HELP HERE
    }
    bool testOutRMS = false;
    if (testOutRMS) {
        if (_simParams->iterationEnd >= 100 && iteration % (_simParams->iterationEnd / 1000) == 0) {
            if (mxRMSETest || myRMSETest || mzRMSETest) {
                std::cout << "Iter. #" << iteration << " | RK" << rkStageName << " ";
                if (mxRMSETest)
                    std::cout << "| RMSE. mx: " << rmse_mx << " ";
                if (myRMSETest)
                    std::cout << "| RMSE. my: " << rmse_my << " ";
                if (mzRMSETest)
                    std::cout << "| RMSE. mz: " << rmse_mz << " ";
                std::cout << std::endl;
            }
            //std::cout << "Iteration #" << iteration <<" | RMSE. mx: " << rmse_mx << " | my: " << rmse_my << " | mz:  " << rmse_mz << std::endl;  // Keep for debugging
        }
    }
    /*
    std::cout << "HERE IN DEMAGFIELD1D: X" << std::endl;
    PrintVector(outDemagX, false);

    std::cout << "HERE IN DEMAGFIELD1D: Y" << std::endl;
    PrintVector(outDemagY, false);

    std::cout << "HERE IN DEMAGFIELD1D: Z" << std::endl;
    PrintVector(outDemagZ, false);
     */

    // EASY FIND
    // Clean-up. Probably could free memory for planDFTForwardMx (etc) earlier in function, but it's safer to be here
    fftw_destroy_plan(planDFTForwardMx);
    fftw_destroy_plan(planDFTForwardMy);
    fftw_destroy_plan(planDFTForwardMz);
    fftw_destroy_plan(planDFTForwardHdX);
    fftw_destroy_plan(planDFTForwardHdY);
    fftw_destroy_plan(planDFTForwardHdZ);
    fftw_destroy_plan(planDFTBackMx);
    fftw_destroy_plan(planDFTBackMy);
    fftw_destroy_plan(planDFTBackMz);
    fftw_destroy_plan(planDFTBackHdX);
    fftw_destroy_plan(planDFTBackHdY);
    fftw_destroy_plan(planDFTBackHdZ);
    fftw_free(mX);
    fftw_free(mY);
    fftw_free(mZ);
    fftw_free(hdX);
    fftw_free(hdY);
    fftw_free(hdZ);
}

void DemagnetisationFields::DemagField1DRealGreen(std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ,
                                           std::vector<double>& inMxTerms, std::vector<double>& inMyTerms, std::vector<double>& inMzTerms,
                                           int iteration, std::string rkStageName) {

    int gotNumSpins = GV.GetNumSpins();  // Vectors and arrays defined in function don't include pinned end terms; 4000 sites in length
    int trueNumSpins = GV.GetNumSpins() + 2;  // Vectors and arrays defined out with function include pinned end terms; 4002 sites in length
    const double imagTerm = 0.0, normalizationFactor = 1.0 / static_cast<double>(gotNumSpins);

    auto fftw_alloc_and_check = [](const char* var_name, const int& size) -> double* {
        /*
         * Temporary lambda function for debugging. Keep within `DemagField1D` until debugging complete.
         * Allocates memory for FFTW-suitable arrays.Currently being cautious, so throw exception if allocation fails,
         * and else set memory to zero to ensure data integrity.
         */
        auto *mTerm = static_cast<double*>(fftw_alloc_real(size));
        if (!mTerm) {
            throw std::runtime_error(std::string("Failed to allocate memory for ") + var_name);
        } else {
            std::memset(mTerm, 0, sizeof(double) * size);
            return mTerm;
        }
    };

    auto PrintFFTWVector = [](int numSpins, fftw_complex* vecToPrint, const char* vecName, double testValReal, double testValImag) {
        /*
         * Temporary lambda function for debugging. Keep within `DemagField1D` until debugging complete.
         * Prints the contents of an FFTW-suitable array. Helps find abnormal and unexpected values (compared to my
         * analytical solutions.
         */
        for (int i = 0; i < numSpins; i++) {
            if (vecToPrint[i][0] != testValReal)
                std::cout << vecName << "[" << i << "][0] = " << vecToPrint[i][0] << std::endl;
            else if (vecToPrint[i][1] != testValImag)
                std::cout << vecName << "[" << i << "][1] = " << vecToPrint[i][1] << std::endl;
        }

    };

    // Assign memory for FFTW-suitable arrays for magnetic moments; used during computation and for RMSE calculation
    auto *mX = fftw_alloc_and_check("mx", trueNumSpins);
    auto *mY = fftw_alloc_and_check("my", trueNumSpins);
    auto *mZ = fftw_alloc_and_check("mz", trueNumSpins);

    // Additional memory for demagnetisation field components; used for output; initialise at zero (empty)
    auto *hdX = fftw_alloc_and_check("hdX", trueNumSpins);
    auto *hdY = fftw_alloc_and_check("hdY", trueNumSpins);
    auto *hdZ = fftw_alloc_and_check("hdZ", trueNumSpins);

    auto *Gxx = fftw_alloc_and_check("Gxx", trueNumSpins);
    auto *Gyy = fftw_alloc_and_check("Gyy", trueNumSpins);
    auto *Gzz = fftw_alloc_and_check("Gzz", trueNumSpins);

    // Populate mX, mY, mZ with input terms
    for (int currentSite = 1; currentSite <= gotNumSpins; ++currentSite) {
        mX[currentSite] = inMxTerms[currentSite];
        mY[currentSite] = inMyTerms[currentSite];
        mZ[currentSite] = inMzTerms[currentSite];
    }

    // Create FFTW plans for forward FFT
    fftw_plan planDFTForwardMx = fftw_plan_r2r_1d(trueNumSpins, mX, mX, FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan planDFTForwardMy = fftw_plan_r2r_1d(trueNumSpins, mY, mY, FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan planDFTForwardMz = fftw_plan_r2r_1d(trueNumSpins, mZ, mZ, FFTW_R2HC, FFTW_ESTIMATE);

    // Execute forward FFT
    fftw_execute(planDFTForwardMx);
    fftw_execute(planDFTForwardMy);
    fftw_execute(planDFTForwardMz);

    // Compute Green's function in Fourier space
    for (int k = 0; k < trueNumSpins; ++k) {
        double k_norm = k * 2.0 * M_PI / gotNumSpins;
        Gxx[k] = _simStates->demagFactors.Nx * std::cos(k_norm);
        Gyy[k] = _simStates->demagFactors.Ny * std::cos(k_norm);
        Gzz[k] = _simStates->demagFactors.Nz * std::cos(k_norm);
    }

    // Multiply in Fourier space (H_d = -N * M)
    for (int k = 0; k < trueNumSpins; ++k) {
        hdX[k] = -Gxx[k] * mX[k];
        hdY[k] = -Gyy[k] * mY[k];
        hdZ[k] = -Gzz[k] * mZ[k];
    }

    // Create FFTW plans for inverse FFT
    fftw_plan planDFTBackHdX = fftw_plan_r2r_1d(trueNumSpins, hdX, hdX, FFTW_HC2R, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdY = fftw_plan_r2r_1d(trueNumSpins, hdY, hdY, FFTW_HC2R, FFTW_ESTIMATE);
    fftw_plan planDFTBackHdZ = fftw_plan_r2r_1d(trueNumSpins, hdZ, hdZ, FFTW_HC2R, FFTW_ESTIMATE);

    // Execute inverse FFT
    fftw_execute(planDFTBackHdX);
    fftw_execute(planDFTBackHdY);
    fftw_execute(planDFTBackHdZ);

    // Normalize results
    for (int currentSite = 1; currentSite <= gotNumSpins; ++currentSite) {
        hdX[currentSite] *= normalizationFactor;
        hdY[currentSite] *= normalizationFactor;
        hdZ[currentSite] *= normalizationFactor;
    }

    // Compute Root Mean Square Error (RMSE) for each component of the magnetic moment
    double rmse_mx = 0, rmse_my = 0, rmse_mz = 0;
    for (int currentSite = 0; currentSite < gotNumSpins; currentSite++) {
        int scaled_len = currentSite + 1;  // mx.size() is 4000 whereas mxTerms.size() is  4002
        // Finding SUM_{i=1}^{N}(data_in[i] - data_reconstructed[i])**2
        rmse_mx += pow(inMxTerms[scaled_len] - mX[currentSite], 2);  // CHECKED: rmse_mx is equal to zero
        rmse_my += pow(inMyTerms[scaled_len] - mY[currentSite], 2);  // CHECKED: rmse_my is equal to zero
        rmse_mz += pow(inMzTerms[scaled_len] - mZ[currentSite], 2);  // CHECKED: rmse_mz is equal to zero
    }

    // Finishing RMSE calculation: RMSE = Sqrt(SUM_{i=1}^{N} * 1/N)
    rmse_mx = sqrt(rmse_mx / gotNumSpins);  // CHECKED: rmse_mx is equal to zero
    rmse_my = sqrt(rmse_my / gotNumSpins);  // CHECKED: rmse_my is equal to zero
    rmse_mz = sqrt(rmse_mz / gotNumSpins);  // CHECKED: rmse_mz is equal to zero


    // FFT error estimation based on machine epsilon and FFTW's scale factor
    const double machineEpsilon = std::numeric_limits<double>::epsilon();
    // const double fftError = machineEpsilon * gotNumSpins * sqrt(gotNumSpins);
    // std::cout << fftError << std::endl;std::exit(0);

    bool mxRMSETest, myRMSETest, mzRMSETest;
    if (fabs(rmse_mx) < machineEpsilon) { mxRMSETest = true; }
    else {mxRMSETest = false;}
    if (fabs(rmse_my) < machineEpsilon) { myRMSETest = true; }
    else { myRMSETest = false; }
    if (fabs(rmse_mz) < machineEpsilon) { mzRMSETest = true; }
    else { mzRMSETest = false; }

    double combinedRMSE = sqrt((pow(rmse_mx, 2) + pow(rmse_my, 2) + pow(rmse_mz, 2)) / 3.0 );
    bool applyDemag;
    if (combinedRMSE < machineEpsilon) {applyDemag = true;}
    else {applyDemag = false;}
    // if (mxRMSETest || myRMSETest || mzRMSETest) {applyDemag = true;}
    // else {applyDemag = false;}

    for (int trueSite = 0; trueSite < gotNumSpins + 2; trueSite++) {
        int scale_site = trueSite - 1;
        if (trueSite == 0 || trueSite == (trueNumSpins - 1) ) {
            // Boundary conditions of output arrays must always be zero; do this to ensure data integrity
            outDemagX[trueSite] = 0.0;
            outDemagY[trueSite] = 0.0;
            outDemagZ[trueSite] = 0.0;
            continue;
        }
        /*
         * outDemagX/outDemagY/outDemagZ are 4002 in length while mX/mY/mZ are 4000 in length. mX (etc) are unmutated throughout
         * the FFT, which inMxTerms (etc) are the original values. If the absolute difference between inMxTerms and mX (etc)
         * is greater than their rmse_mx (etc), then the overwriting of outDemagX (etc) by hdX (etc) should be stopped. This
         * is because the values are likely to just be noise. This is a very crude way of doing this, but it works for now.
         *
         * Caution! outDemagX (etc) are references so overwriting their elements changes what the main program sees.
         */
        outDemagX[trueSite] = applyDemag ? hdX[scale_site] : 0.0; // hdX[scale_site][0]; //(fabs((inMxTerms[trueSite] - mX[scale_site][0])) > fabs(rmse_mx)) ? hdX[scale_site][0] : 0.0;  // HELP HERE
        outDemagY[trueSite] = applyDemag ? hdY[scale_site] : 0.0;// hdY[scale_site][0]; //(fabs((inMyTerms[trueSite] - mY[scale_site][0])) > fabs(rmse_my)) ? hdY[scale_site][0] : 0.0;  // HELP HERE
        outDemagZ[trueSite] = applyDemag ? hdZ[scale_site] : 0.0;// hdZ[scale_site][0]; //(fabs((inMzTerms[trueSite] - mZ[scale_site][0])) > fabs(rmse_mz)) ? hdZ[scale_site][0] : 0.0;  // HELP HERE
    }
    bool testOutRMS = false;
    if (testOutRMS) {
        if (_simParams->iterationEnd >= 100 && iteration % (_simParams->iterationEnd / 1000) == 0) {
            if (mxRMSETest || myRMSETest || mzRMSETest) {
                std::cout << "Iter. #" << iteration << " | RK" << rkStageName << " ";
                if (mxRMSETest)
                    std::cout << "| RMSE. mx: " << rmse_mx << " ";
                if (myRMSETest)
                    std::cout << "| RMSE. my: " << rmse_my << " ";
                if (mzRMSETest)
                    std::cout << "| RMSE. mz: " << rmse_mz << " ";
                std::cout << std::endl;
            }
            //std::cout << "Iteration #" << iteration <<" | RMSE. mx: " << rmse_mx << " | my: " << rmse_my << " | mz:  " << rmse_mz << std::endl;  // Keep for debugging
        }
    }
    /*
    std::cout << "HERE IN DEMAGFIELD1D: X" << std::endl;
    PrintVector(outDemagX, false);

    std::cout << "HERE IN DEMAGFIELD1D: Y" << std::endl;
    PrintVector(outDemagY, false);

    std::cout << "HERE IN DEMAGFIELD1D: Z" << std::endl;
    PrintVector(outDemagZ, false);
     */

    // EASY FIND
    // Clean-up. Probably could free memory for planDFTForwardMx (etc) earlier in function, but it's safer to be here
    fftw_destroy_plan(planDFTForwardMx);
    fftw_destroy_plan(planDFTForwardMy);
    fftw_destroy_plan(planDFTForwardMz);
    fftw_destroy_plan(planDFTBackHdX);
    fftw_destroy_plan(planDFTBackHdY);
    fftw_destroy_plan(planDFTBackHdZ);
    fftw_free(mX);
    fftw_free(mY);
    fftw_free(mZ);
    fftw_free(hdX);
    fftw_free(hdY);
    fftw_free(hdZ);
    fftw_free(Gxx);
    fftw_free(Gyy);
    fftw_free(Gzz);
}

void DemagnetisationFields::DemagFieldsUsingDipoles(std::vector<double> mxTerms, std::vector<double> myTerms,
                                                  std::vector<double> mzTerms, std::vector<int> sitePositions,
                                                  std::vector<double>& outDemagX, std::vector<double>& outDemagY, std::vector<double>& outDemagZ) {
        // Initialization
    std::vector<double> totalDipoleTerms = {0.0, 0.0, 0.0};
    std::vector<double> totalDemagTerms = {0.0, 0.0, 0.0};

    // Copy over common variables and constants
    double exchangeStiffness = 5.3e-17;
    double dipoleValue;

    for (int i = 0; i < mxTerms.size(); i++) {
        // Current site's magnetic moment
        std::vector<double> originSite = {mxTerms[i], myTerms[i], mzTerms[i]};

        for (int j = 0; j < mxTerms.size(); j++) {
            if (i == j) continue; // Skip self-interactions

            // Displacement vector between sites i and j
            std::vector<int> positionVector = {(sitePositions[j] - sitePositions[i]),
                                                  0,
                                                  0};

            double positionVector_norm = std::sqrt(std::pow(positionVector[0], 2) + std::pow(positionVector[1], 2) + std::pow(positionVector[2], 2));
            double positionVector_cubed = std::pow(positionVector_norm, 3);

            // Magnetic moment of site j
            std::vector<double> influencingSite = {mxTerms[j], myTerms[j], mzTerms[j]};

            // Compute scalar product between originSite and influencingSite
            double dotProduct = originSite[0]*influencingSite[0] + originSite[1]*influencingSite[1] + originSite[2]*influencingSite[2];

            // Compute interaction between originSite and influencingSite
            for (int p = 0; p < 3; p++) {
                dipoleValue = 3 * positionVector[p] * dotProduct / positionVector_cubed - influencingSite[p] / positionVector_cubed;
                totalDipoleTerms[p] += dipoleValue;
            }
        }

        // The demag field for each magnetic moment will be the negative of its dipolar field.
        for (int p = 0; p < 3; p++) {
            totalDemagTerms[p] = -totalDipoleTerms[p];
        }
    }
}

float DemagnetisationFields::_calculateDemagFactorsUniformPrism( const double &a, const double &b, const double &c) {
    // Use double during computation for precision, but output as float (values will be of order
    // 1e-6 <= demagnetisationFactor <= 1

    // Use a small offset value to prevent x/0 and log(x) from occurring
    double epsilon = std::numeric_limits<double>::epsilon();

    double r = a * a + b * b + c * c;
    double sqrt_r = std::sqrt(r);

    // Use the exact formula as given in the image
    double demagnetisationFactor = 0.0;

    double ab = a * b;
    double abc = a * b * c;

    double sqrt_ab = std::sqrt(a * a + b * b);
    double sqrt_bc = std::sqrt(b * b + c * c);
    double sqrt_ac = std::sqrt(a * a + c * c);

    demagnetisationFactor += ((b * b - c * c) / (2.0 * b * c)) * std::log((sqrt_r - a) / (sqrt_r + a + epsilon));
    demagnetisationFactor += ((a * a - c * c) / (2.0 * a * c)) * std::log((sqrt_r - b) / (sqrt_r + b + epsilon));
    demagnetisationFactor += (b / (2.0 * c)) * std::log((sqrt_ab + a) / (sqrt_ab - a + epsilon));
    demagnetisationFactor += (a / (2.0 * c)) * std::log((sqrt_ab + b) / (sqrt_ab - b + epsilon));
    demagnetisationFactor += (c / (2.0 * a)) * std::log((sqrt_bc - b) / (sqrt_bc + b + epsilon));
    demagnetisationFactor += (c / (2.0 * b)) * std::log((sqrt_ac - a) / (sqrt_ac + a + epsilon));
    demagnetisationFactor += 2.0 * std::atan(ab / (c * sqrt_r + epsilon));
    demagnetisationFactor += (a * a * a + b * b * b - 2.0 * c * c * c) / (3.0 * abc);
    demagnetisationFactor += ((a * a + b * b - 2.0 * c * c) / (3.0 * abc)) * sqrt_r;
    demagnetisationFactor += (c / ab) * (sqrt_ac + sqrt_bc);
    demagnetisationFactor -= (std::pow((a * a + b * b), 1.5) + std::pow((b * b + c * c), 1.5) + std::pow((a * a + c * c), 1.5)) / (3.0 * abc);
    demagnetisationFactor /= M_PI;

    return static_cast<float>(demagnetisationFactor);
}

void DemagnetisationFields::calculateOneDimension( const std::vector<double> &mxTermsIn,
                                                   const std::vector<double> &myTermsIn,
                                                   const std::vector<double> &mzTermsIn,
                                                   std::vector<std::atomic<double>> &demagFieldXOut,
                                                   std::vector<std::atomic<double>> &demagFieldYOut,
                                                   std::vector<std::atomic<double>> &demagFieldZOut,
                                                   const bool &shouldUseTBB ) {
    // This function is used for parallel calculations. Useful in large systems or when H_eff is complex

    if ( shouldUseTBB ) {
        // tbb::global_control c( tbb::global_control::max_allowed_parallelism, 1 );
        tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins + 1),
            [&](const tbb::blocked_range<int>& range) {
                for (int site = range.begin(); site < range.end(); site++) {
                    std::array<double, 3> tempDemagResults = {0.0, 0.0, 0.0};

                    // Default formula calculates Nz
                    tempDemagResults[0] = -1 * _simStates->demagFactors.Nx * mxTermsIn[site];// / _simParams->satMag;
                    tempDemagResults[1] = -1 * _simStates->demagFactors.Ny * myTermsIn[site];// / _simParams->satMag;
                    tempDemagResults[2] = -1 * _simStates->demagFactors.Nz * mzTermsIn[site];// / _simParams->satMag;

                    // Only change the x-component of the bias field as this is the only component that will be driven dynamically
                    demagFieldXOut[site].fetch_add(tempDemagResults[0]);
                    demagFieldYOut[site].fetch_add(tempDemagResults[1]);
                    demagFieldZOut[site].fetch_add(tempDemagResults[2]);
                }
        }, tbb::auto_partitioner());

    } else {
        throw std::invalid_argument("calculateOneDimension for exchange fields hasn't got CUDA implementation yet");
    }
}