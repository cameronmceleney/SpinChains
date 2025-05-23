//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/DemagnetisationFields.h"

DemagnetisationFields::DemagnetisationFields(SimulationParameters* sharedSimParams,
                                             SimulationStates* sharedSimStates,
                                             SimulationFlags* sharedSimFlags)

   : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}

// New code after refactoring is below
void DemagnetisationFields::initialise(const bool &useDebug) {
    auto length = static_cast<float>(_simParams->systemTotalSpins * _simParams->latticeConstant * 1e9);
    float width = 2.0f;
    float thickness = 8.0f;

    if (_simStates->demagFactors.isDemagFactorDefaultState()) {
        _simStates->demagFactors.Nz = _calculateDemagFactorsUniformPrism(length, width, thickness);
        _simStates->demagFactors.Ny = _calculateDemagFactorsUniformPrism(thickness, length, width);
        _simStates->demagFactors.Nx = _calculateDemagFactorsUniformPrism(width, thickness, length);
    }

    if (!_simStates->demagFactors.isDemagFactorErrorAcceptable()) {
        throw std::runtime_error("DemagnetisationFields::initialise: Demagnetisation factors are not within acceptable bounds");
    }

    if (useDebug) {
        std::cout << _simStates->demagFactors.Nx << ", " << _simStates->demagFactors.Ny << ", "
                              << _simStates->demagFactors.Nz << std::endl;
    }
}

template void DemagnetisationFields::calculateOneDimension<double>(
    const std::vector<double>& mxTerms, const std::vector<double>& myTerms,
    const std::vector<double>& mzTerms, std::vector<double>& demagFieldXOut,
    std::vector<double>& demagFieldYOut, std::vector<double>& demagFieldZOut,
    const DemagMethod& demagMethod, CommonStructures::Parallelisations parallelFlag);

template void DemagnetisationFields::calculateOneDimension<std::atomic<double>>(
    const std::vector<double>& mxTerms, const std::vector<double>& myTerms,
    const std::vector<double>& mzTerms, std::vector<std::atomic<double>>& demagFieldXOut,
    std::vector<std::atomic<double>>& demagFieldYOut, std::vector<std::atomic<double>>& demagFieldZOut,
    const DemagMethod& demagMethod, CommonStructures::Parallelisations parallelFlag);

template <typename T>
void DemagnetisationFields::_addDemagField(int site, const CommonStructures::Vector3D& demagField,
                                           std::vector<T>& demagFieldXOut, std::vector<T>& demagFieldYOut,
                                           std::vector<T>& demagFieldZOut) {
    demagFieldXOut[site] += demagField[0];
    demagFieldYOut[site] += demagField[1];
    demagFieldZOut[site] += demagField[2];
}

template <>
void DemagnetisationFields::_addDemagField<std::atomic<double>>(int site, const CommonStructures::Vector3D& demagField,
                                           std::vector<std::atomic<double>>& demagFieldXOut,
                                           std::vector<std::atomic<double>>& demagFieldYOut,
                                           std::vector<std::atomic<double>>& demagFieldZOut) {
    demagFieldXOut[site].fetch_add(demagField[0], std::memory_order_relaxed);
    demagFieldYOut[site].fetch_add(demagField[1], std::memory_order_relaxed);
    demagFieldZOut[site].fetch_add(demagField[2], std::memory_order_relaxed);
}

template <typename T>
void DemagnetisationFields::calculateOneDimension(const std::vector<double>& mxTerms, const std::vector<double>& myTerms,
                                                  const std::vector<double>& mzTerms, std::vector<T>& demagFieldXOut,
                                                  std::vector<T>& demagFieldYOut, std::vector<T>& demagFieldZOut,
                                                  const DemagnetisationFields::DemagMethod& demagMethod,
                                                  CommonStructures::Parallelisations parallelFlag) {

    if (parallelFlag == CommonStructures::Parallelisations::Multithreaded) {
        tbb::parallel_for(tbb::blocked_range<int>(0, _simParams->systemTotalSpins),
            [&](const tbb::blocked_range<int>& range) {
                for (int site = range.begin(); site < range.end(); site++) {
                    auto tempFieldLocal = _calculateFieldForSite(site, mxTerms, myTerms, mzTerms, demagMethod);
                    _addDemagField(site, tempFieldLocal, demagFieldXOut, demagFieldYOut, demagFieldZOut);
                }
            }, tbb::auto_partitioner());
    } else if (parallelFlag == CommonStructures::Parallelisations::Sequential) {
        for (int site = 0; site < _simParams->systemTotalSpins; site++) {
            auto tempField = _calculateFieldForSite(site, mxTerms, myTerms, mzTerms, demagMethod);
            _addDemagField(site, tempField, demagFieldXOut, demagFieldYOut, demagFieldZOut);
        }
    } else {
        throw std::invalid_argument("DemagnetisationFields::calculateOneDimension: Invalid parallelisation flag");
    }
}


CommonStructures::Vector3D DemagnetisationFields::_calculateFieldForSite(const int &site, const std::vector<double>& mxTerms,
                                                                    const std::vector<double>& myTerms,
                                                                    const std::vector<double>& mzTerms,
                                                                    const DemagnetisationFields::DemagMethod& demagMethod) {
    if (demagMethod == DemagMethod::Demag) {
        return calculateDemagField(1, site, mxTerms, myTerms, mzTerms);
    } else if (demagMethod == DemagMethod::DipoleGreenFunctionReal) {
        auto [dipoleX, dipoleY, dipoleZ] = calculateDipoleFieldGreenReal(1, site, mxTerms, myTerms, mzTerms);
        return {dipoleX[site], dipoleY[site], dipoleZ[site]};
    } else if (demagMethod == DemagMethod::DipoleGreenFunctionComplex) {
        return calculateDipoleFieldGreenComplex(1, site, mxTerms, myTerms, mzTerms);
    } else if (demagMethod == DemagMethod::DipoleBruteForce) {
        return calculateDipoleFieldBruteForce(1, site, mxTerms, myTerms, mzTerms);
    } else {
        throw std::invalid_argument("DemagnetisationFields::calculateOneDimension::calculateFieldForSite Invalid demagnetisation method");
    }
}

CommonStructures::Vector3D DemagnetisationFields::calculateDemagField(int dimension, int site,
                                                                 const std::vector<double>& mxTerms,
                                                                 const std::vector<double>& myTerms,
                                                                 const std::vector<double>& mzTerms) {
    return _calculateDemagField(site, mxTerms, myTerms, mzTerms);
}

CommonStructures::Vector3D DemagnetisationFields::calculateDipoleFieldBruteForce(int dimension, int site,
                                                                 const std::vector<double>& mxTerms,
                                                                 const std::vector<double>& myTerms,
                                                                 const std::vector<double>& mzTerms) {
    return _calculateDipoleFieldBruteForce(site, mxTerms, myTerms, mzTerms);
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
DemagnetisationFields::calculateDipoleFieldGreenReal( int dimension, int site,
                                                      const std::vector<double>& mxTerms,
                                                      const std::vector<double>& myTerms,
                                                      const std::vector<double>& mzTerms) {
    return _calculateDipoleFieldGreen1DReal(site, mxTerms, myTerms, mzTerms);
}

CommonStructures::Vector3D DemagnetisationFields::calculateDipoleFieldGreenComplex( int dimension, int site,
                                                                            const std::vector<double>& mxTerms,
                                                                            const std::vector<double>& myTerms,
                                                                            const std::vector<double>& mzTerms) {
    return _calculateDipoleFieldGreen1DComplex(site, mxTerms, myTerms, mzTerms);
}

CommonStructures::Vector3D DemagnetisationFields::_calculateDemagField(int site, const std::vector<double>& mxTerms, const std::vector<double>& myTerms, const std::vector<double>& mzTerms) {
    return {
        -1 * _simStates->demagFactors.Nx * mxTerms[site],
        -1 * _simStates->demagFactors.Ny * myTerms[site],
        -1 * _simStates->demagFactors.Nz * mzTerms[site]
    };
}

CommonStructures::Vector3D DemagnetisationFields::_calculateDipoleFieldBruteForce(const int &site, const std::vector<double> &mxTerms, const std::vector<double> &myTerms,
                                                                             const std::vector<double> &mzTerms) {
    // Initialization
    CommonStructures::Vector3D localTotalDipoleField{0.0, 0.0, 0.0};

    std::vector<double> originSite = {mxTerms[site], myTerms[site], mzTerms[site]};

    for (int i = 0; i < mxTerms.size(); i++) {
        if (i == site) continue; // Skip self-interactions

        std::vector<int> positionVector = {(i - site), 0, 0};

        double positionVector_norm = std::sqrt(std::pow(positionVector[0], 2) + std::pow(positionVector[1], 2) + std::pow(positionVector[2], 2));
        double positionVector_cubed = std::pow(positionVector_norm, 3);

        // Magnetic moment of site j
        std::vector<double> influencingSite = {mxTerms[i], myTerms[i], mzTerms[i]};

        // Compute scalar product between originSite and influencingSite
        double dotProduct = originSite[0]*influencingSite[0] + originSite[1]*influencingSite[1] + originSite[2]*influencingSite[2];

        // Compute interaction between originSite and influencingSite
        for (int p = 0; p < 3; p++) {
            localTotalDipoleField[p] += 3 * positionVector[p] * dotProduct / positionVector_cubed - influencingSite[p] / positionVector_cubed;
        }
    }

    // The dipolar field contains a negative sign which we include at this final stage (for each component)
    return {-localTotalDipoleField[0], -localTotalDipoleField[1], -localTotalDipoleField[2]};
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

template <typename T>
T* DemagnetisationFields::fftw_alloc_and_check(const char* var_name, int size) {
    T* mTerm = static_cast<T*>(fftw_malloc(sizeof(T) * size));
    if (!mTerm) {
        throw std::runtime_error(std::string("Failed to allocate memory for ") + var_name);
    }
    std::memset(mTerm, 0, sizeof(T) * size);
    return mTerm;
}

void DemagnetisationFields::populateMemory(fftw_complex* mX, fftw_complex* mY, fftw_complex* mZ, const std::vector<double>& inMxTerms,
                                           const std::vector<double>& inMyTerms, const std::vector<double>& inMzTerms,
                                           bool& shouldSkipBounds) {
    int iMin = shouldSkipBounds ? 1 : 0;
    for (int currentSite = iMin; currentSite < _simParams->systemTotalSpins + shouldSkipBounds; ++currentSite) {
        mX[currentSite][0] = inMxTerms[currentSite];
        mX[currentSite][1] = 0.0;
        mY[currentSite][0] = inMyTerms[currentSite];
        mY[currentSite][1] = 0.0;
        mZ[currentSite][0] = inMzTerms[currentSite];
        mZ[currentSite][1] = 0.0;
    }
}

void DemagnetisationFields::populateMemory(double* mX, double* mY, double* mZ, const std::vector<double>& inMxTerms,
                                           const std::vector<double>& inMyTerms, const std::vector<double>& inMzTerms,
                                           bool& shouldSkipBounds) {
    int iMin = shouldSkipBounds ? 1 : 0;
    for (int currentSite = iMin; currentSite < _simParams->systemTotalSpins + shouldSkipBounds; ++currentSite) {
        mX[currentSite] = inMxTerms[currentSite];
        mY[currentSite] = inMyTerms[currentSite];
        mZ[currentSite] = inMzTerms[currentSite];
    }
}

void
DemagnetisationFields::createAndExecuteFFTWPlans( fftw_complex *planXIn, fftw_complex *planYIn, fftw_complex *planZIn,
                                                  fftw_complex *planXOut, fftw_complex *planYOut,
                                                  fftw_complex *planZOut, int trueNumSpins, bool forward ) {
    fftw_plan planX, planY, planZ;

    if (forward) {
        planX = fftw_plan_dft_1d(trueNumSpins, planXIn, planXOut, FFTW_FORWARD, FFTW_ESTIMATE);
        planY = fftw_plan_dft_1d(trueNumSpins, planYIn, planYOut, FFTW_FORWARD, FFTW_ESTIMATE);
        planZ = fftw_plan_dft_1d(trueNumSpins, planZIn, planZOut, FFTW_FORWARD, FFTW_ESTIMATE);
    } else {
        planX = fftw_plan_dft_1d(trueNumSpins, planXOut, planXIn, FFTW_BACKWARD, FFTW_ESTIMATE);
        planY = fftw_plan_dft_1d(trueNumSpins, planYOut, planYIn, FFTW_BACKWARD, FFTW_ESTIMATE);
        planZ = fftw_plan_dft_1d(trueNumSpins, planZOut, planZIn, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (!planX || !planY || !planZ) {
        throw std::runtime_error("DemagnetisationFields::createAndExecuteFFTWPlans -> Failed to create FFTW plans for demagnetisation field calculation");
    }

    fftw_execute(planX);
    fftw_execute(planY);
    fftw_execute(planZ);
}

void DemagnetisationFields::computeGreenFunctionInFourierSpace( double* Gxx, double* Gyy, double* Gzz, int trueNumSpins) {
    double kx, ky, kz;
    double MU_0 = 4 * M_PI * 1e-7;

    for (int k = 1; k < trueNumSpins; ++k) {
        kx = (k < trueNumSpins / 2) ? (2.0 * M_PI * k / trueNumSpins) : (2.0 * M_PI * (k - trueNumSpins) / trueNumSpins);
        ky = kz = 0.0; // In 1D case, only kx is non-zero

        double k_squared = kx * kx + ky * ky + kz * kz;
        if (k_squared > 1e-10) { // Avoid division by zero
            double factor = (MU_0 / (4 * M_PI)) * (1.0 / k_squared);

            Gxx[k] = factor * (kx * kx);
            Gyy[k] = factor * (ky * ky);
            Gzz[k] = factor * (kz * kz);
        } else {
            Gxx[k] = Gyy[k] = Gzz[k] = 0.0;
        }
    }
}

void DemagnetisationFields::computeGreenFunctionInFourierSpace( fftw_complex * Gxx, fftw_complex* Gyy, fftw_complex* Gzz, int trueNumSpins) {

    for (int site = 1; site < trueNumSpins; ++site) {
        double k_norm = site * 2.0 * M_PI / trueNumSpins;
        Gxx[site][0] = 1.0 - cos(k_norm);
        Gxx[site][1] = -sin(k_norm);
        Gyy[site][0] = 1.0 - cos(k_norm);
        Gyy[site][1] = -sin(k_norm);
        Gzz[site][0] = 1.0 - cos(k_norm);
        Gzz[site][1] = -sin(k_norm);
    }
}

void DemagnetisationFields::multiplyInFourierSpace( double* hdX, double* hdY, double* hdZ, const double* Gxx,
                                                    const double* Gyy, const double* Gzz, const fftw_complex *mX,
                                                    const fftw_complex *mY, const fftw_complex *mZ, int trueNumSpins) {
    for (int site = 1; site < trueNumSpins; ++site) {
        hdX[site] = -Gxx[site] * mX[site][0];
        hdY[site] = -Gyy[site] * mY[site][0];
        hdZ[site] = -Gzz[site] * mZ[site][0];
    }
}

void DemagnetisationFields::multiplyInFourierSpace( fftw_complex * hdX, fftw_complex* hdY, fftw_complex* hdZ, const fftw_complex* Gxx, const fftw_complex* Gyy, const fftw_complex* Gzz, const fftw_complex* mX, const fftw_complex* mY, const fftw_complex* mZ, int trueNumSpins) {
    for (int site = 1; site < trueNumSpins; ++site) {
        hdX[site][0] = -Gxx[site][0] * mX[site][0] + Gxx[site][1] * mX[site][1];
        hdX[site][1] = -Gxx[site][0] * mX[site][1] - Gxx[site][1] * mX[site][0];
        hdY[site][0] = -Gyy[site][0] * mY[site][0] + Gyy[site][1] * mY[site][1];
        hdY[site][1] = -Gyy[site][0] * mY[site][1] - Gyy[site][1] * mY[site][0];
        hdZ[site][0] = -Gzz[site][0] * mZ[site][0] + Gzz[site][1] * mZ[site][1];
        hdZ[site][1] = -Gzz[site][0] * mZ[site][1] - Gzz[site][1] * mZ[site][0];
    }
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
DemagnetisationFields::_calculateDipoleFieldGreen1DReal( const int& site, const std::vector<double>& inMxTerms,
                                                         const std::vector<double>& inMyTerms,
                                                         const std::vector<double>& inMzTerms) {

    int trueNumSpins = _simParams->systemTotalSpins + 2;

    const double normalizationFactor = 1.0 / static_cast<int>(_simParams->systemTotalSpins);
    bool skipBounds = true;

    auto* mXIn = fftw_alloc_and_check<fftw_complex>("mxIn", trueNumSpins);
    auto* mX = fftw_alloc_and_check<fftw_complex>("mx", trueNumSpins);
    auto* mYIn = fftw_alloc_and_check<fftw_complex>("myIn", trueNumSpins);
    auto* mY = fftw_alloc_and_check<fftw_complex>("my", trueNumSpins);
    auto* mZIn = fftw_alloc_and_check<fftw_complex>("mzIn", trueNumSpins);
    auto* mZ = fftw_alloc_and_check<fftw_complex>("mz", trueNumSpins);
    auto* hdX = fftw_alloc_and_check<fftw_complex>("hdX", trueNumSpins);
    auto* hdY = fftw_alloc_and_check<fftw_complex>("hdY", trueNumSpins);
    auto* hdZ = fftw_alloc_and_check<fftw_complex>("hdZ", trueNumSpins);
    auto* Gxx = fftw_alloc_and_check<fftw_complex>("Gxx", trueNumSpins);
    auto* Gyy = fftw_alloc_and_check<fftw_complex>("Gyy", trueNumSpins);
    auto* Gzz = fftw_alloc_and_check<fftw_complex>("Gzz", trueNumSpins);

    std::vector<double> resultX(trueNumSpins, 0.0);
    std::vector<double> resultY(trueNumSpins, 0.0);
    std::vector<double> resultZ(trueNumSpins, 0.0);

    try {

        populateMemory(mXIn, mYIn, mZIn, inMxTerms, inMyTerms, inMzTerms, skipBounds);

        createAndExecuteFFTWPlans(mXIn, mYIn, mZIn, mX, mY, mZ, trueNumSpins, true);

        computeGreenFunctionInFourierSpace(Gxx, Gyy, Gzz, trueNumSpins);

        multiplyInFourierSpace(hdX, hdY, hdZ, Gxx, Gyy, Gzz, mX, mY, mZ, trueNumSpins);

        createAndExecuteFFTWPlans(hdX, hdY, hdZ, hdX, hdY, hdZ, trueNumSpins, false);

        for (int i = 0; i < trueNumSpins; ++i) {
            resultX[i] = hdX[i][0] * normalizationFactor;
            resultY[i] = hdY[i][0] * normalizationFactor;
            resultZ[i] = hdZ[i][0] * normalizationFactor;
        }

    } catch (...) {
        fftw_free(mXIn);
        fftw_free(mYIn);
        fftw_free(mZIn);
        fftw_free(mX);
        fftw_free(mY);
        fftw_free(mZ);
        fftw_free(hdX);
        fftw_free(hdY);
        fftw_free(hdZ);
        fftw_free(Gxx);
        fftw_free(Gyy);
        fftw_free(Gzz);
        throw;
    }

    fftw_free(mXIn);
    fftw_free(mYIn);
    fftw_free(mZIn);
    fftw_free(mX);
    fftw_free(mY);
    fftw_free(mZ);
    fftw_free(hdX);
    fftw_free(hdY);
    fftw_free(hdZ);
    fftw_free(Gxx);
    fftw_free(Gyy);
    fftw_free(Gzz);

    return std::make_tuple(resultX, resultY, resultZ);
}

CommonStructures::Vector3D DemagnetisationFields::_calculateDipoleFieldGreen1DComplex(const int& site, const std::vector<double>& inMxTerms,
                                                                                 const std::vector<double>& inMyTerms,
                                                                                 const std::vector<double>& inMzTerms) {

    int gotNumSpins = GV.GetNumSpins();
    int trueNumSpins = gotNumSpins + 2;
    const double normalizationFactor = 1.0 / gotNumSpins;
    bool skipBounds = true;

    auto* mX = fftw_alloc_and_check<fftw_complex>("mx", gotNumSpins);
    auto* mY = fftw_alloc_and_check<fftw_complex>("my", gotNumSpins);
    auto* mZ = fftw_alloc_and_check<fftw_complex>("mz", gotNumSpins);
    auto* hdX = fftw_alloc_and_check<fftw_complex>("hdX", gotNumSpins);
    auto* hdY = fftw_alloc_and_check<fftw_complex>("hxY", gotNumSpins);
    auto* hdZ = fftw_alloc_and_check<fftw_complex>("hdZ", gotNumSpins);
    auto* Gxx = fftw_alloc_and_check<fftw_complex>("Gxx", trueNumSpins);
    auto* Gyy = fftw_alloc_and_check<fftw_complex>("Gyy", trueNumSpins);
    auto* Gzz = fftw_alloc_and_check<fftw_complex>("Gzz", trueNumSpins);

    CommonStructures::Vector3D result = {0.0, 0.0, 0.0};

    try {
        populateMemory(mX, mY, mZ, inMxTerms, inMyTerms, inMzTerms, skipBounds);

        fftw_plan planDFTForwardMx = fftw_plan_dft_1d(trueNumSpins, mX, mX, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan planDFTForwardMy = fftw_plan_dft_1d(trueNumSpins, mY, mY, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan planDFTForwardMz = fftw_plan_dft_1d(trueNumSpins, mZ, mZ, FFTW_FORWARD, FFTW_ESTIMATE);

        fftw_execute(planDFTForwardMx);
        fftw_execute(planDFTForwardMy);
        fftw_execute(planDFTForwardMz);

        computeGreenFunctionInFourierSpace(Gxx, Gyy, Gzz, trueNumSpins);
        multiplyInFourierSpace(hdX, hdY, hdZ, Gxx, Gyy, Gzz, mX, mY, mZ, trueNumSpins);

        fftw_plan planDFTBackHdX = fftw_plan_dft_1d(trueNumSpins, hdX, hdX, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan planDFTBackHdY = fftw_plan_dft_1d(trueNumSpins, hdY, hdY, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan planDFTBackHdZ = fftw_plan_dft_1d(trueNumSpins, hdZ, hdZ, FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(planDFTBackHdX);
        fftw_execute(planDFTBackHdY);
        fftw_execute(planDFTBackHdZ);

        result[0] = hdX[1][0] * normalizationFactor;
        result[1] = hdY[1][0] * normalizationFactor;
        result[2] = hdZ[1][0] * normalizationFactor;
    } catch (...) {
        fftw_free(mX);
        fftw_free(mY);
        fftw_free(mZ);
        fftw_free(hdX);
        fftw_free(hdY);
        fftw_free(hdZ);
        fftw_free(Gxx);
        fftw_free(Gyy);
        fftw_free(Gzz);
        throw;
    }

    fftw_free(mX);
    fftw_free(mY);
    fftw_free(mZ);
    fftw_free(hdX);
    fftw_free(hdY);
    fftw_free(hdZ);
    fftw_free(Gxx);
    fftw_free(Gyy);
    fftw_free(Gzz);

    return result;
}

