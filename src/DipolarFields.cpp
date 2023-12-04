//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/DipolarFields.h"

DipolarFields::DipolarFields( SimulationParameters *sharedSimParams,
                              SimulationStates *sharedSimStates,
                              SimulationFlags *sharedSimFlags )

        : _simParams(sharedSimParams), _simStates(sharedSimStates), _simFlags(sharedSimFlags) {}


double DipolarFields::dipolarKernel3D( const int &originSite, const int &influencingSite, const double &A,
                                       const double &alpha ) {
    // This function is used to calculate the dipolar interaction between two sites. The kernel is defined as:
    // K = 1 / (4 * pi * r^3) * (3 * cos(theta)^2 - 1)
    // where r is the distance between the two sites, and theta is the angle between the two sites.

    // ################################ Declare initial Values ################################
    double dipoleKernelDirect = 0.0, distanceBetweenSites = 0.0, theta = 0.0, cosTheta = 0.0, sinTheta = 0.0;
    double x1 = 0.0, y1 = 0.0, z1 = 0.0, x2 = 0.0, y2 = 0.0, z2 = 0.0;
    // ################################ Calculate distance between sites ################################
    // Calculate the distance between the two sites
    distanceBetweenSites = std::abs(originSite - influencingSite);

    // ################################ Calculate theta ################################
    // Calculate the angle between the two sites
    x1 = originSite;
    y1 = 0.0;
    z1 = 0.0;
    x2 = influencingSite;
    y2 = 0.0;
    z2 = 0.0;

    theta = acos(
            (x1 * x2 + y1 * y2 + z1 * z2) / (sqrt(x1 * x1 + y1 * y1 + z1 * z1) * sqrt(x2 * x2 + y2 * y2 + z2 * z2)));

    // ################################ Calculate cos(theta) ################################
    // Calculate the cosine of the angle between the two sites
    cosTheta = cos(theta);

    // ################################ Calculate sin(theta) ################################
    // Calculate the sine of the angle between the two sites
    sinTheta = sin(theta);

    // ################################ Calculate dipole kernel ################################
    // Calculate the dipole kernel
    dipoleKernelDirect = 1.0 / (4.0 * M_PI * pow(distanceBetweenSites, 3.0));// * (3.0 * pow(cosTheta, 2.0) - 1.0);
    double dipoleKernelIndirect = A * std::exp(-alpha * distanceBetweenSites);
    double dipoleKernel = dipoleKernelDirect + dipoleKernelIndirect;
    return dipoleKernelDirect;
}

double
DipolarFields::dipolarKernel1D( const int &originSite, const int &influencingSite, const std::string &component ) {
    // ################################ Declare initial Values ################################
    double exchangeStiffness = 5.3e-17;

    // ################################ Calculate distance between sites ################################

    if ( _simStates->exchangeVec[influencingSite] == 0 ) {
        // Guard clause to ensure that the exchange vector is not zero
        return 0.0;
    }

    double latticeConstant = std::sqrt(exchangeStiffness / _simStates->exchangeVec[influencingSite]);

    if ( std::isinf(latticeConstant)) {
        // Guard clause to ensure that the lattice constant is not infinite (backup test / temporary)
        throw std::runtime_error(std::string("Lattice constant is infinite!"));
    }

    double positionVector = (influencingSite - originSite) * latticeConstant;
    // ################################ Calculate dipole kernel ################################
    if ( component == "X" ) { return _simParams->PERM_FREESPACE / (2 * M_PI * pow(positionVector, 3.0)); }
    else { throw std::runtime_error(std::string("Invalid component passed to dipolarKernel1D")); }
}

void DipolarFields::DipolarInteraction1D( std::vector<double> inMxTerms, std::vector<double> &outDipoleX ) {
    // Modelling a one-dimensional chain of spins and calculating the dipolar interaction among them.
    // ################################ Declare initial Values ################################
    int trueNumSpins = GV.GetNumSpins() +
                       2;  // Vectors and arrays defined out with function include pinned end terms; 4002 sites in length instead of 4000
    const double imagTerm = 0.0, normalizationFactor = 1.0 / static_cast<double>(GV.GetNumSpins());

    for ( int i = 0; i < trueNumSpins; i++ ) {
        inMxTerms[i] *= _simParams->PERMEABILITY_IRON;
    }

    // ################################ Lambda Functions ################################
    auto fftw_alloc_and_check = []( const char *var_name, const int &size ) -> fftw_complex * {
        /*
         * Temporary lambda function for debugging. Keep within `DemagField1D` until debugging complete.
         * Allocates memory for FFTW-suitable arrays.Currently being cautious, so throw exception if allocation fails,
         * and else set memory to zero to ensure data integrity.
         */
        auto *mTerm = static_cast<fftw_complex *>(fftw_alloc_complex(size));
        if ( !mTerm ) {
            throw std::runtime_error(std::string("Failed to allocate memory for ") + var_name);
        } else {
            std::memset(mTerm, 0, sizeof(fftw_complex) * size);
            return mTerm;
        }
    };

    // ################################ Begin FFT ################################
    // Assign memory for FFTW-suitable arrays for magnetic moments; used during computation and for RMSE calculation
    auto *mX = fftw_alloc_and_check("mX", trueNumSpins);

    auto *mXTransformed = fftw_alloc_and_check("mXTransformed", trueNumSpins);

    // Additional memory for demagnetisation field components; used for output; initialise at zero (empty)
    auto *dipolarKernelX = fftw_alloc_and_check("dipolarKernelX", trueNumSpins);

    auto *HDipoleX = fftw_alloc_and_check("HDipoleX", trueNumSpins);

    // Population of memory for FFTW-suitable arrays
    for ( int currentSite = 0; currentSite < trueNumSpins; currentSite++ ) {
        if ( currentSite == 0 || currentSite == (trueNumSpins - 1)) {
            // Boundary conditions. Hereafter, skip boundary sites during loops
            mX[currentSite][0] = 0.0;
            mX[currentSite][1] = 0.0;
            continue;
        }
        mX[currentSite][0] = inMxTerms[currentSite];
        mX[currentSite][1] = 0.0;
    }

    // Calculate dipolar interaction kernel
    for ( int currentSite = 1; currentSite < (trueNumSpins - 1); currentSite++ ) {
        double sumKernelX = 0.0;
        for ( int influencingSite = 1; influencingSite < (trueNumSpins - 1); influencingSite++ ) {
            if ( currentSite != influencingSite ) {  // Exclude self-interaction
                sumKernelX += dipolarKernel1D(currentSite, influencingSite, "X");
            }
        }
        dipolarKernelX[currentSite][0] = sumKernelX;
        dipolarKernelX[currentSite][1] = 0.0;
    }

    // ################################ Fourier-Space Transformation ################################
    // Create FFTW plans **after** population of memory
    fftw_plan planDFTForwardMx = fftw_plan_dft_1d(trueNumSpins, mX, mXTransformed, FFTW_FORWARD, FFTW_ESTIMATE);
    if ( !planDFTForwardMx ) { throw std::runtime_error("Failed to create planDFTForwardMx"); }

    fftw_plan planDFTForwardHDipoleX = fftw_plan_dft_1d(trueNumSpins, dipolarKernelX, HDipoleX, FFTW_FORWARD,
                                                        FFTW_ESTIMATE);
    if ( !planDFTForwardHDipoleX ) { throw std::runtime_error("Failed to create planDFTForwardHDipoleX"); }

    // Execute the plans **after** their definitions
    fftw_execute(planDFTForwardMx);
    fftw_execute(planDFTForwardHDipoleX);

    // Destroy the plans **after** their execution to save memory in runtime
    fftw_destroy_plan(planDFTForwardMx);
    fftw_destroy_plan(planDFTForwardHDipoleX);

    // ################################ Calculation in Fourier Space ################################
    for ( int currentSite = 1; currentSite < trueNumSpins - 1; currentSite++ ) {
        double tempRealX = mXTransformed[currentSite][0] * HDipoleX[currentSite][0] -
                           mXTransformed[currentSite][1] * HDipoleX[currentSite][1];

        double tempImagX = mXTransformed[currentSite][0] * HDipoleX[currentSite][1] +
                           mXTransformed[currentSite][1] * HDipoleX[currentSite][0];

        mXTransformed[currentSite][0] = tempRealX;
        mXTransformed[currentSite][1] = tempImagX;
    }
    // ################################ Inverse Fourier-Space Transformation ################################
    // Create inverse FFTW plans **after** population of memory
    fftw_plan planDFTBackHDipoleX = fftw_plan_dft_1d(trueNumSpins, mXTransformed, HDipoleX, FFTW_BACKWARD,
                                                     FFTW_ESTIMATE);
    if ( !planDFTBackHDipoleX ) { throw std::runtime_error("Failed to create planDFTBackHDipoleX"); }
    // Execute the plans **after** their definitions
    fftw_execute(planDFTBackHDipoleX);

    // Destroy the plans **after** their execution to save memory in runtime
    fftw_destroy_plan(planDFTBackHDipoleX);

    for ( int currentSite = 0; currentSite < trueNumSpins; currentSite++ ) {
        if ( currentSite == 0 || currentSite == (trueNumSpins - 1)) {
            // Boundary conditions of output arrays must always be zero; do this to ensure data integrity
            outDipoleX[currentSite] = 0.0;
            continue;
        }
        outDipoleX[currentSite] = HDipoleX[currentSite][0] * normalizationFactor;
    }

    // Clean-up. Probably could free memory for planDFTForwardMx (etc) earlier in function, but it's safer to be here
    fftw_free(mX);
    fftw_free(mXTransformed);
    fftw_free(dipolarKernelX);
    fftw_free(HDipoleX);
}

std::vector<double> DipolarFields::DipolarInteractionClassic( std::vector<double> mxTerms, std::vector<double> myTerms,
                                                              std::vector<double> mzTerms,
                                                              std::vector<int> sitePositions ) {

    std::vector<double> totalDipoleTerms = {0.0, 0.0, 0.0};  // Returns Dipole terms for a single site
    // sitePositions contains 3 elements: {site to the left, current site, site to the right}

    for ( int i = 0; i < mxTerms.size(); i++ ) {
        mxTerms[i] *= _simParams->PERMEABILITY_IRON;
        myTerms[i] *= _simParams->PERMEABILITY_IRON;
        mzTerms[i] *= _simParams->PERMEABILITY_IRON;
    }

    double exchangeStiffness = 5.3e-17;
    double exchangeValue;

    // Reference moments
    std::vector<double> originSite = {mxTerms[1], myTerms[1], mzTerms[1]};

    for ( int i = 0; i < mxTerms.size(); i++ ) {
        if ( i == 1 ) {
            // Guard clause to ensure that the origin site is not included in the calculation
            continue;
        }

        if ( _simStates->exchangeVec[sitePositions[i]] == 0 ) {
            // Guard clause to ensure that the exchange vector is not zero
            continue;
        }

        double latticeConstant = std::sqrt(exchangeStiffness / _simStates->exchangeVec[sitePositions[i]]);

        if ( std::isinf(latticeConstant)) {
            // Guard clause to ensure that the lattice constant is not infinite (backup test / temporary)
            continue;
        }

        std::vector<double> positionVector = {(sitePositions[i] - sitePositions[1]) * latticeConstant, 0, 0};

        double positionVector_norm = std::sqrt(std::pow(positionVector[0], 2) + std::pow(positionVector[1], 2)
                                               + std::pow(positionVector[2], 2));

        double positionVector_cubed = std::pow(positionVector_norm, 3);
        double positionVector_fifth = std::pow(positionVector_norm, 5);

        // Moment at site i
        std::vector<double> influencingSite = {mxTerms[i], myTerms[i], mzTerms[i]};

        // double originSiteDotInfluencingSite = originSite[0] * influencingSite[0] + originSite[1] * influencingSite[1] + originSite[2] * influencingSite[2];
        double originSiteDotPosition = originSite[0] * positionVector[0] + originSite[1] * positionVector[1] +
                                       originSite[2] * positionVector[2];
        double influencingSiteDotPosition =
                influencingSite[0] * positionVector[0] + influencingSite[1] * positionVector[1] +
                influencingSite[2] * positionVector[2];

        for ( int j = 0; j < 3; j++ ) {
            double DipoleValue = _simParams->dipoleConstant *
                                 ((3.0 * positionVector[j] * influencingSiteDotPosition) / positionVector_fifth -
                                  influencingSite[j] / positionVector_cubed);
            totalDipoleTerms[j] += DipoleValue;
        }
    }
    return totalDipoleTerms;
}

std::vector<double> DipolarFields::DipolarInteractionIntralayer( std::vector<std::vector<double>> &mTerms,
                                                                 int &currentSite, const int &currentLayer,
                                                                 const double &exchangeStiffness ) {
    /* This function calculates the dipolar interaction between the current site and its neighbours within a single layer.
     *
     * WARNING. This function assumes that every site is aligned along the x-axis which is only valid for specific
     * spin chains. This function will need to be modified to account for arbitrary spin chains.
     *
     */
    std::vector<double> totalDipoleTerms = {0.0, 0.0, 0.0};

    int vecLength, originIndex;
    if ( _simParams->numNeighbours == 0 ) {
        // Guard clause to ensure that the number of neighbours is not zero
        return totalDipoleTerms;
    } else if ( _simParams->numNeighbours < 0 ) {
        vecLength = _simStates->layerSpinsInChain[currentLayer];
        originIndex = currentSite - _simParams->numSpinsInABC - 1;
    } else {
        vecLength = 2 * _simParams->numNeighbours + 1;
        originIndex = vecLength / 2 + 1;
    }

    if ( vecLength < 0 )
        std::cout << "Error: vecLength is less than zero" << std::endl;

    // Could combine these to be a single vector for memory improvements
    std::vector<double> mxTerms(vecLength, 0);
    std::vector<double> myTerms(vecLength, 0);
    std::vector<double> mzTerms(vecLength, 0);
    std::vector<int> sitePositions(vecLength, 0);

    /* This IF statement will be optimised away by passing an array of the form [x1,x2,...,y1,y2,...,z1,z2,...] when
     * CUDA is implemented; instead of giving a general 2D mTerms vector and then forcing this function to flatten.
     */
    int iFV = 0; // index flat vector
    if ( _simParams->numNeighbours < 0 ) {
        for ( int site = _simParams->numSpinsInABC + 1; site <= vecLength + _simParams->numSpinsInABC; site++ ) {
            // Flatting the vectors
            mxTerms[iFV] = mTerms[site][0] * _simParams->PERMEABILITY_IRON;
            myTerms[iFV] = mTerms[site][1] * _simParams->PERMEABILITY_IRON;
            mzTerms[iFV] = mTerms[site][2] * _simParams->PERMEABILITY_IRON;
            sitePositions[iFV] = site;
            iFV++;
        }
    } else {
        for ( int site = currentSite - _simParams->numNeighbours;
              site <= currentSite + _simParams->numNeighbours; site++ ) {
            if ( site < _simParams->numSpinsInABC or
                 site >= _simStates->layerSpinsInChain[currentLayer] + _simParams->numSpinsInABC ) {
                // Guard clause to skip trying assignment of any element when the index is negative
                continue;
            }
            // Flatting the vectors
            mxTerms[iFV] = mTerms[site][0] * _simParams->PERMEABILITY_IRON;
            myTerms[iFV] = mTerms[site][1] * _simParams->PERMEABILITY_IRON;
            mzTerms[iFV] = mTerms[site][2] * _simParams->PERMEABILITY_IRON;
            sitePositions[iFV] = site;
            iFV++;
        }
    }
    // Here to improve readability; could be removed to improve performance
    std::vector<double> originSite = {mxTerms[originIndex], myTerms[originIndex], mzTerms[originIndex]};

    // Start of the loop over the neighbours
    for ( int i = 0; i < vecLength; i++ ) {
        if ( i == originIndex ) {
            // Guard clause to ensure that the origin site is not included in the calculation
            continue;
        }

        // Moment at site i. Here to improve readability; could be removed to improve performance
        std::vector<double> influencingSite = {mxTerms[i], myTerms[i], mzTerms[i]};
        if ( influencingSite[0] == 0.0 && influencingSite[1] == 0.0 && influencingSite[2] == 0.0 ) {
            // If influencing site components are all zero, then they don't impact the calculation. So can be skipped
            continue;
        }

        if ( exchangeStiffness == 0.0 || _simStates->exchangeVec[sitePositions[i] - 1] == 0.0 ) {
            // simParams->exchangeVec[sitePositions[i]-1] refers to exchange vector to the LHS of the current site; [i] is RHS
            continue;
        }

        double latticeConstant = std::sqrt(exchangeStiffness / _simStates->exchangeVec[sitePositions[i] - 1]);

        if ( std::isinf(latticeConstant)) {
            // Guard clause to ensure that the lattice constant is not infinite (backup test / temporary)
            continue;
        }

        std::vector<double> positionVector = {(sitePositions[i] - sitePositions[originIndex]) * latticeConstant, 0, 0};

        double positionVector_norm = positionVector[0];  // Simplifies to this for only a single component

        double positionVector_cubed = std::pow(positionVector_norm, 3);
        double positionVector_fifth = std::pow(positionVector_norm, 5);

        if ( positionVector_cubed == 0.0 || positionVector_fifth == 0.0 ) {
            // Could use an epsilon value here to avoid division by zero and to make the code more efficient
            continue;
        }
        // Calculate the dot products
        double originSiteDotPosition = originSite[0] * positionVector[0];

        double influencingSiteDotPosition = influencingSite[0] * positionVector[0];

        for ( int j = 0; j < 3; j++ ) {
            // Calculate the dipole-dipole coupling term
            double DipoleValue = _simParams->dipoleConstant * (((3.0 * positionVector[j] * influencingSiteDotPosition)
                                                                / positionVector_fifth) -
                                                               influencingSite[j] / positionVector_cubed);
            totalDipoleTerms[j] += DipoleValue;
        }
    }

    return totalDipoleTerms;
}

std::vector<double> DipolarFields::DipolarInteractionInterlayer( std::vector<std::vector<double>> &mTermsLayer1,
                                                                 std::vector<std::vector<double>> &mTermsLayer2,
                                                                 int &currentSite, const int &currentLayer,
                                                                 const int &otherLayer ) {
    std::vector<double> totalDipoleTerms = {0.0, 0.0, 0.0};
    bool findAdj = false;

    double exchangeStiffness = 5.3e-17;
    double interlayerExchange = 132.0;  // Interlayer exchange coupling in Tesla

    if ( currentSite <= _simParams->numSpinsInABC or
         currentSite > (_simStates->layerSpinsInChain[currentLayer] + _simParams->numSpinsInABC)) {
        return {0.0, 0.0, 0.0};  // Ensure currentSite is valid within the current (target) layer
    }

    // Calculate the dipolar coupling for chain1
    std::vector<double> totalDipoleTermsLayer1 = DipolarInteractionIntralayer(mTermsLayer1, currentSite, currentLayer,
                                                                              exchangeStiffness);

    std::vector<double> totalDipoleTermsOtherChains;
    if ( findAdj ) {
        totalDipoleTermsOtherChains = DipolarInteractionInterlayerAdjacent(mTermsLayer1, mTermsLayer2,
                                                                           _simParams->numNeighbours, currentSite,
                                                                           currentLayer, exchangeStiffness,
                                                                           interlayerExchange);
    } else {
        totalDipoleTermsOtherChains = DipolarInteractionInterlayerAll(mTermsLayer1, mTermsLayer2,
                                                                      currentSite, currentLayer, otherLayer,
                                                                      exchangeStiffness, interlayerExchange);
    }

    // Finally add the three dipole terms to get the total dipole term for a site in chain 1
    for ( int i = 0; i < 3; i++ ) {
        totalDipoleTerms[i] += totalDipoleTermsLayer1[i] + totalDipoleTermsOtherChains[i];
    }

    return totalDipoleTerms;
}

std::vector<double> DipolarFields::DipolarInteractionInterlayerAll( std::vector<std::vector<double>> &mTermsLayer1,
                                                                    std::vector<std::vector<double>> &mTermsLayer2,
                                                                    int &currentSite, const int &currentLayer,
                                                                    const int &otherLayer, double &exchangeStiffness,
                                                                    double &interlayerExchange ) {
    /* Calculate the dipolar interaction between a site in Layer1 (chain 1), and every other site in another layer (chain 2).
     *
     * WARNING. This function is only valid for the following conditions: the two layers are parallel; the distance
     * between sites in each layer is the same; there is no z-component involved in the position coordinates. The
     * removal of the z-coordinate allows for fewer calculations.
     *
     */

    std::vector<double> totalDipolarInteractionInterlayer = {0.0, 0.0, 0.0};

    // Stop-gap code to prevent memory-access violation error. Needs fixed in the future
    int chainTwoOffset;
    if ( !_simFlags->shouldDriveAllLayers ) {
        chainTwoOffset = _simStates->layerSpinsInChain[otherLayer] + _simParams->numSpinsInABC;
    } else { chainTwoOffset = _simStates->layerSpinsInChain[currentLayer] + _simParams->numSpinsInABC; }

    for ( int otherSite = 0; otherSite < mTermsLayer2.size(); otherSite++ ) {
        if ( otherSite > _simParams->numSpinsInABC and otherSite <= chainTwoOffset ) {
            // Exclude damped regions as they are aphysical and will lead to incorrect results

            double intralayerLatticeConstant = std::sqrt(exchangeStiffness / _simStates->exchangeVec[currentSite]);
            double interlayerLatticeConstant = std::sqrt(exchangeStiffness / interlayerExchange);

            if ( std::isinf(intralayerLatticeConstant) or std::isinf(interlayerLatticeConstant)) {
                // Guard clause to ensure that the lattice constant is not infinite (backup test / temporary)
                continue;
            }

            std::vector<double> positionVector = {(otherSite - currentSite) * intralayerLatticeConstant,
                                                  interlayerLatticeConstant, 0};

            double positionVector_norm = std::sqrt(std::pow(positionVector[0], 2) + std::pow(positionVector[1], 2));
            double positionVector_cubed = std::pow(positionVector_norm, 3);
            double positionVector_fifth = std::pow(positionVector_norm, 5);

            std::vector<double> originSite = {mTermsLayer1[currentSite][0] * _simParams->PERMEABILITY_IRON,
                                              mTermsLayer1[currentSite][1] * _simParams->PERMEABILITY_IRON,
                                              mTermsLayer1[currentSite][2] * _simParams->PERMEABILITY_IRON};
            std::vector<double> influencingSite = {mTermsLayer2[otherSite][0] * _simParams->PERMEABILITY_IRON,
                                                   mTermsLayer2[otherSite][1] * _simParams->PERMEABILITY_IRON,
                                                   mTermsLayer2[otherSite][2] * _simParams->PERMEABILITY_IRON};

            double originSiteDotPosition = originSite[0] * positionVector[0] + originSite[1] * positionVector[1];
            double influencingSiteDotPosition = influencingSite[0] * positionVector[0]
                                                + influencingSite[1] * positionVector[1];

            for ( int j = 0; j < 3; j++ ) {
                double DipoleValue =
                        _simParams->dipoleConstant * (((3.0 * positionVector[j] * influencingSiteDotPosition)
                                                       / positionVector_fifth) -
                                                      influencingSite[j] / positionVector_cubed);
                totalDipolarInteractionInterlayer[j] += DipoleValue;
            }

        }
    }

    return totalDipolarInteractionInterlayer;
}

std::vector<double> DipolarFields::DipolarInteractionInterlayerAdjacent( std::vector<std::vector<double>> &mTermsChain1,
                                                                         std::vector<std::vector<double>> &mTermsChain2,
                                                                         int &numNeighbours, int &currentSite,
                                                                         const int &currentLayer,
                                                                         double &exchangeStiffness,
                                                                         double &interlayerExchange ) {
    /* Calculate the dipolar interaction between a site in Layer1 (chain 1), and every other site in another layer
     * (chain 2) within the driving region.
     *
     * WARNING. This function is only valid for the following conditions: the two layers are parallel; the distance
     * between sites in each layer is the same; there is no x- or z-components involved in the position coordinates; the driving
     * region of Layer1 overlaps exactly with the intended dipolar driven region of Layer2.
     *
     * To calculate the dipolar interaction between every site in another chain and your current site, use
     * `DipolarInteractionInterlayerAll`
     *
     */

    std::vector<double> totalDipolarInteractionInterlayer = {0.0, 0.0, 0.0};

    // Stop-gap code to prevent memory-access violation error. Needs fixed in the future
    int chainTwoOffset;
    if ( !_simFlags->shouldDriveAllLayers ) {
        chainTwoOffset = _simStates->layerSpinsInChain[0] + _simParams->numSpinsInABC;
    } else { chainTwoOffset = _simStates->layerSpinsInChain[currentLayer] + _simParams->numSpinsInABC; }

    // Check if currentSite is a valid index for mTermsChain2 before calculations
    if ( currentSite > _simParams->numSpinsInABC and currentSite <= chainTwoOffset ) {
        // Could also calculate coupling for each site in chain 2, but this is computationally expensive

        double interlayerLatticeConstant = std::sqrt(exchangeStiffness / interlayerExchange);

        std::vector<double> positionVector = {0, interlayerLatticeConstant, 0};
        double positionVector_norm = positionVector[1];  // Simplifies to this for only a single component
        double positionVector_cubed = std::pow(positionVector_norm, 3);
        double positionVector_fifth = std::pow(positionVector_norm, 5);

        std::vector<double> originSite = {mTermsChain1[currentSite][0] * _simParams->PERMEABILITY_IRON,
                                          mTermsChain1[currentSite][1] * _simParams->PERMEABILITY_IRON,
                                          mTermsChain1[currentSite][2] * _simParams->PERMEABILITY_IRON};

        std::vector<double> influencingSite = {mTermsChain2[currentSite][0] * _simParams->PERMEABILITY_IRON,
                                               mTermsChain2[currentSite][1] * _simParams->PERMEABILITY_IRON,
                                               mTermsChain2[currentSite][2] * _simParams->PERMEABILITY_IRON};

        double originSiteDotPosition = originSite[1] * positionVector[1];
        double influencingSiteDotPosition = influencingSite[1] * positionVector[1];

        for ( int j = 0; j < 3; j++ ) {
            double DipoleValue = _simParams->dipoleConstant *
                                 ((3.0 * positionVector[j] * influencingSiteDotPosition) / positionVector_fifth
                                  - influencingSite[j] / positionVector_cubed);
            totalDipolarInteractionInterlayer[j] += DipoleValue;
        }
    }

    return totalDipolarInteractionInterlayer;
}

std::vector<double>
DipolarFields::DipolarInteractionIntralayerDebug( std::vector<std::vector<double>> &mTerms, int &numNeighbours,
                                                  int &currentSite, const int &currentLayer ) {
    std::vector<double> totalDipoleTerms = {0.0, 0.0, 0.0};

    double exchangeStiffness = 5.3e-17;

    if ( _simFlags->debugFunc ) { std::cout << "DB2.1 | "; }

    int vecLength, originIndex;
    if ( numNeighbours == 0 ) {
        // Guard clause to ensure that the number of neighbours is not zero
        return totalDipoleTerms;
    } else if ( numNeighbours < 0 ) {
        vecLength = _simStates->layerSpinsInChain[currentLayer];
        originIndex = currentSite - _simParams->numSpinsInABC - 1;
    } else {
        vecLength = 2 * numNeighbours + 1;
        originIndex = vecLength / 2 + 1;
    }
    if ( _simFlags->debugFunc ) { std::cout << "DB2.2 | "; }

    if ( vecLength < 0 )
        std::cout << "Error: vecLength is less than zero" << std::endl;

    // Could combine these to be a single vector for memory improvements
    std::vector<double> mxTerms(vecLength, 0);
    std::vector<double> myTerms(vecLength, 0);
    std::vector<double> mzTerms(vecLength, 0);
    std::vector<int> sitePositions(vecLength, 0);
    if ( _simFlags->debugFunc ) { std::cout << "DB2.3 | "; }

    int iFV = 0; // index flat vector
    if ( numNeighbours < 0 ) {
        for ( int site = _simParams->numSpinsInABC + 1; site <= vecLength + _simParams->numSpinsInABC; site++ ) {
            // Flatting the vectors
            mxTerms[iFV] = mTerms[site][0] * _simParams->PERMEABILITY_IRON;
            myTerms[iFV] = mTerms[site][1] * _simParams->PERMEABILITY_IRON;
            mzTerms[iFV] = mTerms[site][2] * _simParams->PERMEABILITY_IRON;
            sitePositions[iFV] = site;
            iFV++;
        }
    } else {
        for ( int site = currentSite - numNeighbours; site <= currentSite + numNeighbours; site++ ) {
            if ( site < _simParams->numSpinsInABC or
                 site >= _simStates->layerSpinsInChain[currentLayer] + _simParams->numSpinsInABC ) {
                // Guard clause to skip trying assignment of any element when the index is negative
                continue;
            }
            // Flatting the vectors
            mxTerms[iFV] = mTerms[site][0] * _simParams->PERMEABILITY_IRON;
            myTerms[iFV] = mTerms[site][1] * _simParams->PERMEABILITY_IRON;
            mzTerms[iFV] = mTerms[site][2] * _simParams->PERMEABILITY_IRON;
            sitePositions[iFV] = site;
            iFV++;
        }
    }
    if ( _simFlags->debugFunc ) { std::cout << "DB2.4 | "; }
    // Here to improve readability; could be removed to improve performance
    std::vector<double> originSite = {mxTerms[originIndex], myTerms[originIndex], mzTerms[originIndex]};
    if ( _simFlags->debugFunc ) {
        std::cout << " INTRA Origin: [" << originSite[0] << ", " << originSite[1] << ", " << originSite[2] << "] | ";
    }

    if ( _simFlags->debugFunc ) { std::cout << "DB2.5 | "; }
    // Start of the loop over the neighbours
    for ( int i = 0; i < vecLength; i++ ) {
        if ( _simFlags->debugFunc ) { std::cout << "\nDB2.5.0 (" << i + 1 << ") | "; }
        if ( i == originIndex ) {
            // Guard clause to ensure that the origin site is not included in the calculation
            if ( _simFlags->debugFunc ) { std::cout << "DB2.5.0 - Skip 1 (Same Site) "; }
            continue;
        }

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.1 | "; }
        // Moment at site i. Here to improve readability; could be removed to improve performance
        std::vector<double> influencingSite = {mxTerms[i], myTerms[i], mzTerms[i]};
        if ( _simFlags->debugFunc ) {
            std::cout << "INTRA Influe: [" << influencingSite[0] << ", " << influencingSite[1] << ", "
                      << influencingSite[2] << "] | ";
        }
        if ( influencingSite[0] == 0.0 && influencingSite[1] == 0.0 && influencingSite[2] == 0.0 ) {
            // If influencing site components are all zero, then they don't impact the calculation. So can be skipped
            if ( _simFlags->debugFunc ) { std::cout << "DB2.5.1 - Skip 2 (All influencing components zero)"; }
            continue;
        }

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.2 | "; }
        if ( exchangeStiffness == 0.0 || _simStates->exchangeVec[sitePositions[i] - 1] == 0.0 ) {
            if ( _simFlags->debugFunc ) { std::cout << "DB2.5.2 - Skip 1 (exchange vector zero)"; }
            continue;
        }

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.3 | "; }
        double latticeConstant = std::sqrt(exchangeStiffness / _simStates->exchangeVec[sitePositions[i] - 1]);

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.4 | "; }
        if ( std::isinf(latticeConstant)) {
            // Guard clause to ensure that the lattice constant is not infinite (backup test / temporary)
            if ( _simFlags->debugFunc ) { std::cout << "DB2.5.4 - Skip 1 (lattice constant inf)"; }
            continue;
        }

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.5 | "; }
        std::vector<double> positionVector = {(sitePositions[i] - sitePositions[originIndex]) * latticeConstant, 0, 0};

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.6 | "; }
        double positionVector_norm = std::sqrt(std::pow(positionVector[0], 2) + std::pow(positionVector[1], 2)
                                               + std::pow(positionVector[2], 2));

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.7 | "; }
        double positionVector_cubed = std::pow(positionVector_norm, 3);
        double positionVector_fifth = std::pow(positionVector_norm, 5);

        if ( positionVector_cubed == 0.0 || positionVector_fifth == 0.0 ) {
            // Could use an epsilon value here to avoid division by zero and to make code more efficient
            if ( _simFlags->debugFunc ) { std::cout << "DB2.5.7 - Skip 1 (position vector ^3 &&/|| ^5 zero)"; }
            continue;
        }
        // Calculate the dot products
        double originSiteDotPosition = originSite[0] * positionVector[0] + originSite[1] * positionVector[1] +
                                       originSite[2] * positionVector[2];

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.8 | "; }
        double influencingSiteDotPosition =
                influencingSite[0] * positionVector[0] + influencingSite[1] * positionVector[1] +
                influencingSite[2] * positionVector[2];

        if ( _simFlags->debugFunc ) { std::cout << "DB2.5.9 | "; }
        for ( int j = 0; j < 3; j++ ) {
            // Calculate the dipole-dipole coupling term
            double DipoleValue = _simParams->dipoleConstant *
                                 ((3.0 * positionVector[j] * influencingSiteDotPosition) / positionVector_fifth -
                                  influencingSite[j] / positionVector_cubed);
            totalDipoleTerms[j] += DipoleValue;
        }
    }
    if ( _simFlags->debugFunc ) { std::cout << "\nDB2.6 | "; }

    return totalDipoleTerms;
}

std::vector<double> DipolarFields::DipolarInteractionInterlayerDebug( std::vector<std::vector<double>> &mTermsChain1,
                                                                      std::vector<std::vector<double>> &mTermsChain2,
                                                                      int &numNeighbours, int &currentSite,
                                                                      const int &currentLayer ) {
    std::vector<double> totalDipoleTerms = {0.0, 0.0, 0.0};

    double exchangeStiffness = 5.3e-17;
    double interlayerExchange = 132.0;  // Interlayer exchange coupling in Tesla

    if ( _simFlags->debugFunc ) { std::cout << "DB1.1 | "; }

    if ( currentSite <= _simParams->numSpinsInABC or
         currentSite > (_simStates->layerSpinsInChain[currentLayer] + _simParams->numSpinsInABC)) {
        return {0.0, 0.0, 0.0};  // Ensure currentSite is a valid index for mTermsChain1
    }
    if ( _simFlags->debugFunc ) { std::cout << "DB1.2 | "; }
    // Stop-gap code to prevent memory-access violation error. Needs fixed in the future
    int testLength;
    if ( !_simFlags->shouldDriveAllLayers ) {
        testLength = _simStates->layerSpinsInChain[0] + _simParams->numSpinsInABC;
    } else { testLength = _simStates->layerSpinsInChain[currentLayer] + _simParams->numSpinsInABC; }
    if ( _simFlags->debugFunc ) { std::cout << "DB1.3 / DB2.0 | "; }
    // Calculate the dipolar coupling for chain1
    std::vector<double> totalDipoleTermsChain1 = DipolarInteractionIntralayer(mTermsChain1, currentSite, currentLayer,
                                                                              exchangeStiffness);
    if ( _simFlags->debugFunc ) { std::cout << "DB1.4 | "; }
    // Check if currentSite is a valid index for mTermsChain2 before calculations
    if ( currentSite > _simParams->numSpinsInABC and currentSite <= testLength ) {
        // Could also calculate coupling for each site in chain 2, but this is computationally expensive

        // Here we use the same calculations as in the original function but for two spins at the same site in different chains
        // Assuming the chains are parallel and the distance between them is latticeConstant
        double interlayerLatticeConstant = std::sqrt(exchangeStiffness / interlayerExchange);
        std::vector<double> positionVector = {0, interlayerLatticeConstant, 0};

        if ( _simFlags->debugFunc ) { std::cout << "DB1.5 | "; }
        double positionVector_norm = std::sqrt(std::pow(positionVector[0], 2) + std::pow(positionVector[1], 2)
                                               + std::pow(positionVector[2], 2));
        double positionVector_cubed = std::pow(positionVector_norm, 3);
        double positionVector_fifth = std::pow(positionVector_norm, 5);

        if ( _simFlags->debugFunc ) { std::cout << "DB1.6 | "; }
        std::vector<double> originSite = {mTermsChain1[currentSite][0] * _simParams->PERMEABILITY_IRON,
                                          mTermsChain1[currentSite][1] * _simParams->PERMEABILITY_IRON,
                                          mTermsChain1[currentSite][2] * _simParams->PERMEABILITY_IRON};
        if ( _simFlags->debugFunc ) {
            std::cout << "INTER Origin: [" << originSite[0] << ", " << originSite[1] << ", " << originSite[2] << "] | ";
        }

        if ( _simFlags->debugFunc ) { std::cout << "DB1.7 | "; }
        std::vector<double> influencingSite = {mTermsChain2[currentSite][0] * _simParams->PERMEABILITY_IRON,
                                               mTermsChain2[currentSite][1] * _simParams->PERMEABILITY_IRON,
                                               mTermsChain2[currentSite][2] * _simParams->PERMEABILITY_IRON};
        if ( _simFlags->debugFunc ) {
            std::cout << "INTER Influe: [" << influencingSite[0] << ", " << influencingSite[1] << ", "
                      << influencingSite[2] << "] | ";
        }

        if ( _simFlags->debugFunc ) { std::cout << "DB1.8 | "; }
        double originSiteDotPosition = originSite[0] * positionVector[0] + originSite[1] * positionVector[1]
                                       + originSite[2] * positionVector[2];
        double influencingSiteDotPosition =
                influencingSite[0] * positionVector[0] + influencingSite[1] * positionVector[1]
                + influencingSite[2] * positionVector[2];

        if ( _simFlags->debugFunc ) { std::cout << "DB1.9 | "; }
        for ( int j = 0; j < 3; j++ ) {
            double DipoleValue = _simParams->dipoleConstant *
                                 ((3.0 * positionVector[j] * influencingSiteDotPosition) / positionVector_fifth
                                  - influencingSite[j] / positionVector_cubed);
            // Only contains y terms so can skip j != 1 (x and z)
            totalDipoleTerms[j] += DipoleValue;
        }
        if ( _simFlags->debugFunc ) { std::cout << "DB1.10 | "; }
    }

    // Finally add the three dipole terms to get the total dipole term for a site in chain 1
    for ( int i = 0; i < 3; i++ ) {
        // Only contains x positions so can skip i > 1 (y & z)
        totalDipoleTerms[i] += totalDipoleTermsChain1[i];
    }
    if ( _simFlags->debugFunc ) { std::cout << "DB1.11 | "; }
    if ( _simFlags->debugFunc && currentLayer == 1 ) {
        std::cout << "totalDipoleTerms: [" << totalDipoleTerms[0] << " " << totalDipoleTerms[1] << " "
                  << totalDipoleTerms[2] << "]";
    }
    return totalDipoleTerms;
}

void DipolarFields::DipolarInteractionClassicThreaded( const int &currentSite, const std::vector<double> &mxTerms,
                                                       const std::vector<double> &myTerms,
                                                       const std::vector<double> &mzTerms,
                                                       double &dipoleXOut, double &dipoleYOut, double &dipoleZOut ) {
    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This should be called within a FOR loop as it only returns single values for the dipole terms

    if ( currentSite <= _simParams->numSpinsInABC ||
         currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
        // Guard clause to ensure that the current site is valid
        dipoleXOut = 0.0;
        dipoleYOut = 0.0;
        dipoleZOut = 0.0;
        return;
    }
    std::array<double, 3> calculationResults = _calculateClassicThreaded(currentSite, mxTerms, myTerms, mzTerms);
    dipoleXOut = calculationResults[0];
    dipoleYOut = calculationResults[1];
    dipoleZOut = calculationResults[2];
}

void DipolarFields::DipolarInteractionClassicThreaded( const std::vector<double> &mxTerms,
                                                       const std::vector<double> &myTerms,
                                                       const std::vector<double> &mzTerms,
                                                       std::vector<double> &dipoleXOut, std::vector<double> &dipoleYOut,
                                                       std::vector<double> &dipoleZOut ) {
    // Note that this function can only be used for a 1D spinchain where the signal is along the x-axis
    // This should be called OUTSIDE a for loop as this method has a nested for/reduce structure

    tbb::parallel_for(tbb::blocked_range<int>(1, _simParams->systemTotalSpins),
                      [&]( const tbb::blocked_range<int> tbbRange ) {
                          for ( int currentSite = tbbRange.begin(); currentSite <= tbbRange.end(); currentSite++ ) {
                              if ( currentSite <= _simParams->numSpinsInABC ||
                                   currentSite > (_simParams->systemTotalSpins + _simParams->numSpinsInABC)) {
                                  // Guard clause to ensure that the current site is valid
                                  dipoleXOut[currentSite] = 0.0;
                                  dipoleYOut[currentSite] = 0.0;
                                  dipoleZOut[currentSite] = 0.0;
                                  return;
                              }

                              std::array<double, 3> calculationResults = _calculateClassicThreaded(currentSite, mxTerms,
                                                                                                   myTerms, mzTerms);

                              dipoleXOut[currentSite] = calculationResults[0];
                              dipoleYOut[currentSite] = calculationResults[1];
                              dipoleZOut[currentSite] = calculationResults[2];
                          }
                      }, tbb::auto_partitioner()
    );
}

std::array<double, 3>
DipolarFields::_calculateClassicThreaded( const int &currentSite, const std::vector<double> &mxTerms,
                                          const std::vector<double> &myTerms,
                                          const std::vector<double> &mzTerms ) {

    struct dipoleSumsLocal {
        double dipoleXAtSite;
        double dipoleYAtSite;
        double dipoleZAtSite;

        dipoleSumsLocal() : dipoleXAtSite(0.0), dipoleYAtSite(0.0), dipoleZAtSite(0.0) {}

        dipoleSumsLocal( double x, double y, double z ) : dipoleXAtSite(x), dipoleYAtSite(y), dipoleZAtSite(z) {}
    };

    dipoleSumsLocal localSums = tbb::parallel_reduce(
            tbb::blocked_range<int>(1, _simParams->systemTotalSpins),  // Intentional loop bounds
            dipoleSumsLocal(),
            [&]( const tbb::blocked_range<int> tbbRange, dipoleSumsLocal localInFor ) -> dipoleSumsLocal {
                for ( int iSite = tbbRange.begin(); iSite <= tbbRange.end(); iSite++ ) {
                    if ( iSite <= _simParams->numSpinsInABC ||
                         iSite > (_simParams->numSpinsInChain + _simParams->numSpinsInABC))
                        continue;  // Ensure that the absorbing regions are not included in the calculation

                    if ( iSite == currentSite )
                        continue;  // Ensure that the origin site isn't included in calculation

                    if ( _simStates->exchangeVec[iSite] == 0.0 )
                        continue;  // Ensure that the exchange vector is not zero

                    double latticeConstant = std::sqrt(_simParams->exchangeStiffness / _simStates->exchangeVec[iSite]);

                    if ( std::isinf(latticeConstant))
                        continue;  // Ensure that the lattice constant is not infinite (backup test / temporary)

                    std::array<double, 3> positionVectorLocal = {(iSite - currentSite) * latticeConstant, 0, 0};

                    double positionNormLocal = std::sqrt(std::pow(positionVectorLocal[0], 2));
                    double positionNormCubed = std::pow(positionNormLocal, 3);
                    double positionNormFifth = std::pow(positionNormLocal, 5);

                    std::array<double, 3> influencingSiteLocal = {mxTerms[iSite], myTerms[iSite], mzTerms[iSite]};

                    // double originSiteDotInfluencingSite = originSite[0] * influencingSite[0];
                    // double originSiteDotPosition = originSite[0] * positionVector[0];
                    double influencingSiteDotPositionLocal = influencingSiteLocal[0] * positionVectorLocal[0];

                    localInFor.dipoleXAtSite += _simParams->dipoleConstant *
                                                ((3.0 * positionVectorLocal[0] * influencingSiteDotPositionLocal) /
                                                 positionNormFifth
                                                 - influencingSiteLocal[0] / positionNormCubed);
                    localInFor.dipoleYAtSite = 0.0; // Intentionally zero for debugging
                    localInFor.dipoleZAtSite = 0.0; // Intentionally zero for debugging
                }
                return localInFor;
            },
            []( const dipoleSumsLocal &a, const dipoleSumsLocal &b ) -> dipoleSumsLocal {
                return {a.dipoleXAtSite + b.dipoleXAtSite, a.dipoleYAtSite + b.dipoleYAtSite,
                        a.dipoleZAtSite + b.dipoleZAtSite};
            }
    );

    return {localSums.dipoleXAtSite, localSums.dipoleYAtSite, localSums.dipoleZAtSite};
}