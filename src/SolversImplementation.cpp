//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/SolversImplementation.h"

SolversImplementation::SolversImplementation( std::shared_ptr<SimulationParameters> sharedSimParams,
                                              std::shared_ptr<SimulationStates> sharedSimStates,
                                              std::shared_ptr<SimulationFlags> sharedSimFlags )

        : SolversSuperClass(std::move(sharedSimParams), std::move(sharedSimStates), std::move(sharedSimFlags)),
          demagField(simParams.get(), simStates.get(), simFlags.get()),
          dmInteraction(simParams.get(), simStates.get(), simFlags.get()),
          effectiveField(simParams.get(), simStates.get(), simFlags.get()),
          dipolarField(simParams.get(), simStates.get(), simFlags.get()),
          llg(simParams.get(), simStates.get(), simFlags.get()) {}

void SolversImplementation::_testShockwaveConditions( double iteration ) {

    if ( simFlags->shouldDriveCeaseEarly ) {
        // and (simParams->isShockwaveOn and simParams->isShockwaveAtMax)) {
        if ( simFlags->isShockwaveOn and not simFlags->isShockwaveAtMax )
            std::cout << "Shock not at maximum when cut-off" << std::endl;

        if ( iteration >= simParams->iterationEnd * simParams->iterEndShock ) {
            // Shockwave begins once simulation is a certain % complete
            simFlags->hasShockwave = false;
            simFlags->isShockwaveOn = false;
            simParams->dynamicBiasField = 0;
        }

        return;
    }

    // If method is triggered, then the applied biasFieldDriving is increased by the scale factor shockwaveScaling
    if ( simFlags->hasShockwave and not simFlags->isShockwaveOn ) {
        if ( iteration >= simParams->iterationEnd * simParams->iterStartShock ) {
            // Shockwave begins once simulation is a certain % complete
            simFlags->isShockwaveOn = true;
            simParams->dynamicBiasField = simParams->shockwaveInitialStrength;
        }

        return;
    }

    if ( simFlags->isShockwaveOn and not simFlags->isShockwaveAtMax ) {
        simParams->dynamicBiasField += simParams->shockwaveStepsize;

        if ( simParams->dynamicBiasField >= simParams->shockwaveMax ) {
            simParams->dynamicBiasField = simParams->shockwaveMax;
            simFlags->isShockwaveAtMax = true;

        }
        return;

    }

}

void SolversImplementation::SolveRK2Classic() {
    std::shared_ptr<SolversDataHandling> childNMData = std::make_shared<SolversDataHandling>(simParams, simStates,
                                                                                             simFlags);
    // Only uses a single spin chain to solve the RK2 midpoint method.

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");

    if ( simFlags->isFerromagnetic ) {
        childNMData->InformUserOfCodeType("RK2 Midpoint (Classic)(FM)");
        childNMData->CreateFileHeader(mxRK2File, "RK2 Midpoint (Classic)(FM)");
    } else if ( !simFlags->isFerromagnetic ) {
        childNMData->InformUserOfCodeType("RK2 Midpoint (Classic)(AFM)");
        childNMData->CreateFileHeader(mxRK2File, "RK2 Midpoint (Classic)(AFM)");
    }

    if ( GV.GetEmailWhenCompleted()) {
        childNMData->CreateMetadata();
    }

    progressbar bar(100);

    std::vector<double> demagX(GV.GetNumSpins() + 2, 0.0), demagY(GV.GetNumSpins() + 2, 0.0), demagZ(
            GV.GetNumSpins() + 2, 0.0);
    //std::vector<double> dipoleX(GV.GetNumSpins() + 2, 0.0), dipoleY(GV.GetNumSpins() + 2, 0.0), dipoleZ(GV.GetNumSpins() + 2, 0.0);

    for ( int iteration = simParams->iterationStart; iteration <= simParams->iterationEnd; iteration++ ) {

        if ( simParams->iterationEnd >= 100 && iteration % (simParams->iterationEnd / 100) == 0 )
            bar.update(); // Doesn't work for fewer than 100 iterations

        _testShockwaveConditions(iteration);

        double t0 = simParams->totalTime;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = simParams->mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(GV.GetNumSpins() + 2, 0), my1(GV.GetNumSpins() + 2, 0), mz1(GV.GetNumSpins() + 2, 0);
        // EASY FIND
        // if (simFlags->hasDipolar)
        //     dipolarField.DipolarInteraction1D(simStates->mx0, dipoleX);
        if ( simFlags->hasDemagIntense ) {
            demagField.DemagnetisationFieldIntense(demagX, demagY, demagZ, simStates->mx0, simStates->my0,
                                                   simStates->mz0);
        } else if ( simFlags->hasDemagFFT ) {
            std::string rkStageName = "2-1";
            demagField.DemagField1DReal(demagX, demagY, demagZ, simStates->mx0, simStates->my0, simStates->mz0,
                                        iteration, rkStageName);

            //std::cout << "Iteration #" << iteration <<" | RMSE. mx: " << rmse_mx << " | my: " << rmse_my << " | mz:  " << rmse_mz << std::endl;  // Keep for debugging


            //if (iteration > 0) {std::cout << "Stage 1" << std::endl; PrintVector(demagZ, false);}
            /*
            std::cout << "HERE IN FUNCTION: X" << std::endl;
            PrintVector(demagX, false);
            std::cout << "HERE IN FUNCTION: Y" << std::endl;
            PrintVector(demagY, false);
            std::cout << "HERE IN FUNCTION: Z" << std::endl;
            PrintVector(demagZ, false);

            if (demagX == demagX2)
                std::cout << "HERE IN demagX are the same" << std::endl;
            else {
                std::cout << "HERE IN demagX are NOT the same" << std::endl;
                for (int i = 0; i < demagX.size(); i++) {
                    if (demagX[i] != demagX2[i])
                        std::cout << i << " | demagX (real)" << demagX[i] << " demagX2 (copy) " << demagX2[i] << std::endl;
                }
            }
            if (demagY == demagY2)
                std::cout << "HERE IN demagY are the same" << std::endl;
            else {
                std::cout << "HERE IN demagY are NOT the same" << std::endl;
                for (int i = 0; i < demagY.size(); i++) {
                    if (demagY[i] != demagY2[i])
                        std::cout << i << " | demagY (real)" << demagY[i] << " demagY2 (copy) " << demagY2[i] << std::endl;
                }
            }
            if (demagZ == demagZ2)
                std::cout << "HERE IN demagZ are the same" << std::endl;
            else {
                std::cout << "HERE IN demagZ are NOT the same" << std::endl;
                for (int i = 0; i < demagZ.size(); i++) {
                    if (demagZ[i] != demagZ2[i])
                        std::cout << i << " | demagZ (real)" << demagZ[i] << " demagZ2 (copy) " << demagZ2[i] << std::endl;
                }
            }
            std::exit(0);
             */
        }

        // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)
        // RK2 Stage 1. Takes initial conditions as inputs.

        for ( int site = 1; site <= GV.GetNumSpins(); site++ ) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;
            /*
            if ((simParams->iterationEnd >= 100 && iteration % (simParams->iterationEnd / 1000) == 0) && site == 500) {
                std::cout << "Iter. #" << iteration << " ";
                std::cout << "| mx: " << simParams->mx0[200] << " - H_dx: " << demagX[200] << " ";
                std::cout << "| my: " << simParams->my0[200] << " - H_dy: " << demagY[200] << " ";
                std::cout << "| mz: " << simParams->mz0[200] << " - H_dz: " << demagZ[200] << " ";
                std::cout << "| mTot: " << sqrt(pow(simParams->mx0[site], 2) + pow(simParams->my0[site], 2) + pow(simParams->mz0[site], 2)) << " ";
                std::cout << std::endl;
                //std::cout << "| mz: " << simParams->mz0[200] << " ";
                //std::cout << "| H_dz: " << demagZ[200] << " ";
                //std::cout << "| mTot: " << sqrt(pow(simParams->mx0[site], 2) + pow(simParams->my0[site], 2) + pow(simParams->mz0[site], 2)) << " ";
                //std::cout << std::endl;
            }
            */
            double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if ( simFlags->hasDipolar ) {
                std::vector<double> mxTermsForDipole = {simStates->mx0[spinLHS], simStates->mx0[site],
                                                        simStates->mx0[spinRHS]};
                std::vector<double> myTermsForDipole = {simStates->my0[spinLHS], simStates->my0[site],
                                                        simStates->my0[spinRHS]};
                std::vector<double> mzTermsForDipole = {simStates->mz0[spinLHS], simStates->mz0[site],
                                                        simStates->mz0[spinRHS]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = dipolarField.DipolarInteractionClassic(mxTermsForDipole,
                                                                                         myTermsForDipole,
                                                                                         mzTermsForDipole,
                                                                                         siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            }

            double dmiZ = 0;
            if ( simFlags->hasDMI ) {
                std::vector<double> mxTermsForDMI = {simStates->mx0[spinLHS], simStates->mx0[site]};
                std::vector<double> myTermsForDMI = {simStates->my0[spinLHS], simStates->my0[site]};
                std::vector<double> mzTermsForDMI = {simStates->mz0[spinLHS], simStates->mz0[site]};
                std::vector<double> dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                              myTermsForDMI,
                                                                              mzTermsForDMI);

                dmiZ = dmiTerms[2];
            }

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK0 = effectiveField.EffectiveFieldX(site, 0, simStates->mx0[spinLHS], simStates->mx0[site],
                                                         simStates->mx0[spinRHS], dipoleX, demagX[site],
                                                         dmiZ, t0);
            double hyK0 = effectiveField.EffectiveFieldY(site, 0, simStates->my0[spinLHS], simStates->my0[site],
                                                         simStates->my0[spinRHS], dipoleY, demagY[site],
                                                         dmiZ);
            double hzK0 = effectiveField.EffectiveFieldZ(site, 0, simStates->mz0[spinLHS], simStates->mz0[site],
                                                         simStates->mz0[spinRHS], dipoleZ, demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK1 = llg.MagneticMomentX(site, simStates->mx0[site], simStates->my0[site], simStates->mz0[site],
                                              hxK0, hyK0, hzK0);
            double myK1 = llg.MagneticMomentY(site, simStates->mx0[site], simStates->my0[site], simStates->mz0[site],
                                              hxK0, hyK0, hzK0);
            double mzK1 = llg.MagneticMomentZ(site, simStates->mx0[site], simStates->my0[site], simStates->mz0[site],
                                              hxK0, hyK0, hzK0);

            // Find (m0 + k1/2) for each site, which is used in the next stage.
            mx1[site] = simStates->mx0[site] + simParams->stepsizeHalf * mxK1;
            my1[site] = simStates->my0[site] + simParams->stepsizeHalf * myK1;
            mz1[site] = simStates->mz0[site] + simParams->stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        // EASY FIND
        std::vector<double> mx2(GV.GetNumSpins() + 2, 0), my2(GV.GetNumSpins() + 2, 0), mz2(GV.GetNumSpins() + 2, 0);
        std::fill(demagX.begin(), demagX.end(), 0.0);
        std::fill(demagY.begin(), demagY.end(), 0.0);
        std::fill(demagZ.begin(), demagZ.end(), 0.0);
        //std::fill(dipoleX.begin(), dipoleX.end(), 0.0); std::fill(dipoleY.begin(), dipoleY.end(), 0.0); std::fill(dipoleZ.begin(), dipoleZ.end(), 0.0);

        //if (simFlags->hasDipolar)
        //    dipolarField.DipolarInteraction1D(mx1, dipoleX);
        if ( simFlags->hasDemagIntense ) {
            demagField.DemagnetisationFieldIntense(demagX, demagY, demagZ, mx1, my1, mz1);
        } else if ( simFlags->hasDemagFFT ) {
            std::string rkStageName = "2-2";
            demagField.DemagField1DReal(demagX, demagY, demagZ, mx1, my1, mz1, iteration, rkStageName);
            // if (iteration > 0) {std::cout << "Stage 2" << std::endl; PrintVector(demagZ, false);}
        }

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
        for ( int site = 1; site <= GV.GetNumSpins(); site++ ) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if ( simFlags->hasDipolar ) {
                std::vector<double> mxTermsForDipole = {mx1[spinLHS], mx1[site], mx1[spinRHS]};
                std::vector<double> myTermsForDipole = {my1[spinLHS], my1[site], my1[spinRHS]};
                std::vector<double> mzTermsForDipole = {mz1[spinLHS], mz1[site], mz1[spinRHS]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = dipolarField.DipolarInteractionClassic(mxTermsForDipole,
                                                                                         myTermsForDipole,
                                                                                         mzTermsForDipole,
                                                                                         siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            }

            double dmiZ = 0;
            if ( simFlags->hasDMI ) {
                std::vector<double> mxTermsForDMI = {simStates->mx0[spinLHS], simStates->mx0[site]};
                std::vector<double> myTermsForDMI = {simStates->my0[spinLHS], simStates->my0[site]};
                std::vector<double> mzTermsForDMI = {simStates->mz0[spinLHS], simStates->mz0[site]};
                std::vector<double> dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                              myTermsForDMI,
                                                                              mzTermsForDMI);

                dmiZ = dmiTerms[2];
            }

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK1 = effectiveField.EffectiveFieldX(site, 0, mx1[spinLHS], mx1[site], mx1[spinRHS], dipoleX,
                                                         demagX[site], dmiZ, t0);
            double hyK1 = effectiveField.EffectiveFieldY(site, 0, my1[spinLHS], my1[site], my1[spinRHS], dipoleY,
                                                         demagY[site], dmiZ);
            double hzK1 = effectiveField.EffectiveFieldZ(site, 0, mz1[spinLHS], mz1[site], mz1[spinRHS], dipoleZ,
                                                         demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK2 = llg.MagneticMomentX(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double myK2 = llg.MagneticMomentY(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double mzK2 = llg.MagneticMomentZ(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);

            mx2[site] = simStates->mx0[site] + simParams->stepsize * mxK2;
            my2[site] = simStates->my0[site] + simParams->stepsize * myK2;
            mz2[site] = simStates->mz0[site] + simParams->stepsize * mzK2;

            if ( simFlags->shouldTrackMagneticMomentNorm ) {
                double mIterationNorm = sqrt(pow(mx2[site], 2) + pow(my2[site], 2) + pow(mz2[site], 2));
                if ((simParams->largestMNorm) > (1.0 - mIterationNorm)) {
                    simParams->largestMNorm = (1.0 - mIterationNorm);
                }
                //if (mIterationNorm > 1.00005) {throw std::runtime_error("mag. moments are no longer below <= 1.00005");}
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        simStates->mx0.clear();
        simStates->my0.clear();
        simStates->mz0.clear();
        mx1.clear();
        my1.clear();
        mz1.clear();

        childNMData->SaveDataToFile(mxRK2File, mx2, iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        simStates->mx0 = mx2;
        simStates->my0 = my2;
        simStates->mz0 = mz2;

        if ( iteration == simParams->forceStopAtIteration )
            exit(0);

        simParams->totalTime += simParams->stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    if ( GV.GetEmailWhenCompleted()) {
        childNMData->CreateMetadata(true);
    }

    if ( simFlags->shouldTrackMagneticMomentNorm )
        std::cout << "\nMax norm. value of M is: " << simParams->largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void SolversImplementation::SolveRK2() {
    std::shared_ptr<SolversDataHandling> childNMData = std::make_shared<SolversDataHandling>(simParams, simStates,
                                                                                             simFlags);
    // Uses multiple layers to solve the RK2 midpoint method. See the documentation for more details.

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    std::ofstream mxRK2File1(GV.GetFilePath() + "rk2_mx1_" + GV.GetFileNameBase() + ".csv");

    // User information and file header is magnetic-material specific.
    if ( simFlags->isFerromagnetic ) {
        childNMData->InformUserOfCodeType("RK2 Midpoint (FM)");
        childNMData->CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)", false, 0);
        childNMData->CreateFileHeader(mxRK2File1, "RK2 Midpoint (FM)", false, 1);
    } else if ( !simFlags->isFerromagnetic ) {
        childNMData->InformUserOfCodeType("RK2 Midpoint (AFM)");
        childNMData->CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
        childNMData->CreateFileHeader(mxRK2File1, "RK2 Midpoint (AFM)");
    }

    if ( GV.GetEmailWhenCompleted()) {
        childNMData->CreateMetadata();
    }

    progressbar bar(100);

    std::vector<double> demagX(GV.GetNumSpins() + 2, 0.0), demagY(GV.GetNumSpins() + 2, 0.0), demagZ(
            GV.GetNumSpins() + 2, 0.0);

    for ( int iteration = simParams->iterationStart; iteration <= simParams->iterationEnd; iteration++ ) {

        if ( simParams->iterationEnd >= 100 && iteration % (simParams->iterationEnd / 100) == 0 )
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        _testShockwaveConditions(iteration);

        double t0 = simParams->totalTime;

        for ( int layer = 0; layer < simParams->numLayers; layer++ ) {
            // RK2 Stage 1. Takes initial conditions as inputs.

            for ( int site = 1; site <= simStates->layerTotalSpins[layer]; site++ ) {
                // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)

                // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
                int spinLHS = site - 1, spinRHS = site + 1;

                double mxLHS = simStates->m0Nest[layer][spinLHS][0], mxMID = simStates->m0Nest[layer][site][0], mxRHS = simStates->m0Nest[layer][spinRHS][0];
                double myLHS = simStates->m0Nest[layer][spinLHS][1], myMID = simStates->m0Nest[layer][site][1], myRHS = simStates->m0Nest[layer][spinRHS][1];
                double mzLHS = simStates->m0Nest[layer][spinLHS][2], mzMID = simStates->m0Nest[layer][site][2], mzRHS = simStates->m0Nest[layer][spinRHS][2];

                double dipoleX, dipoleY, dipoleZ;
                if ( simFlags->hasDipolar ) {

                    int layer1, layer2;
                    if ( layer == 0 ) {
                        layer1 = 0;
                        layer2 = 1;
                    } else if ( layer == 1 ) {
                        layer1 = 1;
                        layer2 = 0;
                    }

                    if ( simFlags->debugFunc ) {
                        std::cout << "\n\niteration: " << iteration << " | layer: " << layer << " | site: " << site
                                  << std::endl;
                    }
                    std::vector<double> dipoleTerms = dipolarField.DipolarInteractionInterlayer(
                            simStates->m0Nest[layer1], simStates->m0Nest[layer2], site,
                            layer1, layer2);

                    dipoleX = dipoleTerms[0];
                    dipoleY = dipoleTerms[1];
                    dipoleZ = dipoleTerms[2];
                } else {
                    dipoleX = 0;
                    dipoleY = 0;
                    dipoleZ = 0;
                }

                double dmiZ = 0;
                if ( simFlags->hasDMI ) {
                    std::vector<double> mxTermsForDMI = {simStates->m0Nest[layer][spinLHS][0],
                                                         simStates->m0Nest[layer][site][0]};
                    std::vector<double> myTermsForDMI = {simStates->m0Nest[layer][spinLHS][1],
                                                         simStates->m0Nest[layer][site][1]};
                    std::vector<double> mzTermsForDMI = {simStates->m0Nest[layer][spinLHS][2],
                                                         simStates->m0Nest[layer][site][2]};
                    std::vector<double> dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                                  myTermsForDMI,
                                                                                  mzTermsForDMI);

                    dmiZ = dmiTerms[2];
                }

                // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                double hxK0 = effectiveField.EffectiveFieldX(site, layer, mxLHS, mxMID, mxRHS, dipoleX, demagX[site],
                                                             dmiZ, t0);
                double hyK0 = effectiveField.EffectiveFieldY(site, layer, myLHS, myMID, myRHS, dipoleY, demagY[site],
                                                             dmiZ);
                double hzK0 = effectiveField.EffectiveFieldZ(site, layer, mzLHS, mzMID, mzRHS, dipoleZ, demagZ[site]);

                // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
                double mxK1 = llg.MagneticMomentX(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
                double myK1 = llg.MagneticMomentY(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
                double mzK1 = llg.MagneticMomentZ(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);

                // Find (m0 + k1/2) for each site, which is used in the next stage.
                simStates->m1Nest[layer][site][0] = mxMID + simParams->stepsizeHalf * mxK1;
                simStates->m1Nest[layer][site][1] = myMID + simParams->stepsizeHalf * myK1;
                simStates->m1Nest[layer][site][2] = mzMID + simParams->stepsizeHalf * mzK1;
            }
        }

        for ( int layer = 0; layer < simParams->numLayers; layer++ ) {
            // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
            for ( int site = 1; site <= simStates->layerTotalSpins[layer]; site++ ) {

                // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
                int spinLHS = site - 1, spinRHS = site + 1;

                double mxLHS = simStates->m1Nest[layer][spinLHS][0], mxMID = simStates->m1Nest[layer][site][0], mxRHS = simStates->m1Nest[layer][spinRHS][0];
                double myLHS = simStates->m1Nest[layer][spinLHS][1], myMID = simStates->m1Nest[layer][site][1], myRHS = simStates->m1Nest[layer][spinRHS][1];
                double mzLHS = simStates->m1Nest[layer][spinLHS][2], mzMID = simStates->m1Nest[layer][site][2], mzRHS = simStates->m1Nest[layer][spinRHS][2];

                double dipoleX, dipoleY, dipoleZ;
                if ( simFlags->hasDipolar ) {

                    int layer1, layer2;
                    if ( layer == 0 ) {
                        layer1 = 0;
                        layer2 = 1;
                    } else if ( layer == 1 ) {
                        layer1 = 1;
                        layer2 = 0;
                    }

                    int debugCounter = 0;  // To make sure debug outputs only occur during the first RK2 stage, not this second stage
                    if ( simFlags->debugFunc ) {
                        simFlags->debugFunc = false;
                        debugCounter++;
                    }
                    std::vector<double> dipoleTerms = dipolarField.DipolarInteractionInterlayer(
                            simStates->m1Nest[layer1], simStates->m1Nest[layer2], site,
                            layer1, layer2);

                    if ( debugCounter > 0 ) { simFlags->debugFunc = true; }

                    dipoleX = dipoleTerms[0];
                    dipoleY = dipoleTerms[1];
                    dipoleZ = dipoleTerms[2];
                } else {
                    dipoleX = 0;
                    dipoleY = 0;
                    dipoleZ = 0;
                }

                double dmiZ = 0;
                if ( simFlags->hasDMI ) {
                    std::vector<double> mxTermsForDMI = {simStates->m0Nest[layer][spinLHS][0],
                                                         simStates->m0Nest[layer][site][0]};
                    std::vector<double> myTermsForDMI = {simStates->m0Nest[layer][spinLHS][1],
                                                         simStates->m0Nest[layer][site][1]};
                    std::vector<double> mzTermsForDMI = {simStates->m0Nest[layer][spinLHS][2],
                                                         simStates->m0Nest[layer][site][2]};
                    std::vector<double> dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                                  myTermsForDMI,
                                                                                  mzTermsForDMI);

                    dmiZ = dmiTerms[2];
                }
                // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                double hxK1 = effectiveField.EffectiveFieldX(site, layer, mxLHS, mxMID, mxRHS, dipoleX, demagX[site],
                                                             dmiZ, t0);
                double hyK1 = effectiveField.EffectiveFieldY(site, layer, myLHS, myMID, myRHS, dipoleY, demagY[site],
                                                             dmiZ);
                double hzK1 = effectiveField.EffectiveFieldZ(site, layer, mzLHS, mzMID, mzRHS, dipoleZ, demagZ[site]);

                // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
                double mxK2 = llg.MagneticMomentX(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
                double myK2 = llg.MagneticMomentY(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
                double mzK2 = llg.MagneticMomentZ(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);

                simStates->m2Nest[layer][site][0] = simStates->m0Nest[layer][site][0] + simParams->stepsize * mxK2;
                simStates->m2Nest[layer][site][1] = simStates->m0Nest[layer][site][1] + simParams->stepsize * myK2;
                simStates->m2Nest[layer][site][2] = simStates->m0Nest[layer][site][2] + simParams->stepsize * mzK2;

                if ( simFlags->shouldTrackMagneticMomentNorm ) {
                    double mIterationNorm = sqrt(
                            pow(simStates->m2Nest[layer][site][0], 2) + pow(simStates->m2Nest[layer][site][1], 2) +
                            pow(simStates->m2Nest[layer][site][2], 2));
                    if ((simStates->largestMNormMulti[layer]) >
                        (1.0 - mIterationNorm)) { simStates->largestMNormMulti[layer] = (1.0 - mIterationNorm); }
                }
            }
        }

        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */

        childNMData->SaveDataToFileMultilayer(mxRK2File, simStates->m2Nest[0], iteration, 0);
        childNMData->SaveDataToFileMultilayer(mxRK2File1, simStates->m2Nest[1], iteration, 1);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        simStates->m0Nest = simStates->m2Nest;

        if ( iteration == simParams->forceStopAtIteration )
            exit(0);

        simParams->totalTime += simParams->stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    if ( GV.GetEmailWhenCompleted()) {
        childNMData->CreateMetadata(true);
    }

    if ( simFlags->shouldTrackMagneticMomentNorm ) {
        std::cout << "\nMax norm. values of M are: ";
        for ( int i = 0; i < simStates->largestMNormMulti.size(); i++ ) {
            std::cout << "Layer " << i << ": " << simStates->largestMNormMulti[i] << " | ";
        }
    }

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void SolversImplementation::RK2Parallel() {
    // Only works for a 1D spin chain
    std::shared_ptr<SolversDataHandling> solverOutputs = std::make_shared<SolversDataHandling>(simParams, simStates,
                                                                                               simFlags);

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxOutputFile(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");

    if ( simFlags->isFerromagnetic ) {
        solverOutputs->InformUserOfCodeType("RK2 Midpoint (Parallel)(FM)");
        solverOutputs->CreateFileHeader(mxOutputFile, "RK2 Midpoint (Parallel)(FM)");
    } else {
        solverOutputs->InformUserOfCodeType("RK2 Midpoint (AFM)");
        solverOutputs->CreateFileHeader(mxOutputFile, "RK2 Midpoint (AFM)");
    }

    if ( GV.GetEmailWhenCompleted()) { solverOutputs->CreateMetadata(); }

    progressbar bar(100);

    // Faster to declare memory and then read/write repeatedly
    // TODO. As this is all multi-threaded, having such larger vectors is excessive. Change methods to work by element (so that onl mx/my/mz are large in memory)
    std::vector<double> demagX(simParams->systemTotalSpins + 2, 0.0), demagY(simParams->systemTotalSpins + 2,
                                                                             0.0), demagZ(
            simParams->systemTotalSpins + 2, 0.0);
    std::vector<double> dipoleX(simParams->systemTotalSpins + 2, 0.0), dipoleY(simParams->systemTotalSpins + 2,
                                                                               0.0), dipoleZ(
            simParams->systemTotalSpins + 2, 0.0);
    std::vector<double> dmiX(simParams->systemTotalSpins + 2, 0.0), dmiY(simParams->systemTotalSpins + 2, 0.0), dmiZ(
            simParams->systemTotalSpins + 2, 0.0);

    for ( int iteration = simParams->iterationStart; iteration <= simParams->iterationEnd; iteration++ ) {

        if ( simParams->iterationEnd >= 100 && iteration % (simParams->iterationEnd / 100) == 0 )
            bar.update(); // Doesn't work for fewer than 100 iterations

        _testShockwaveConditions(iteration);

        double t0 = simParams->totalTime;

        // RK2 Stage 1. Takes initial conditions as inputs. The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = simParams->mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(simParams->systemTotalSpins + 2, 0), my1(simParams->systemTotalSpins + 2, 0), mz1(
                simParams->systemTotalSpins + 2, 0);
        RK2StageMultithreaded(simStates->mx0, simStates->my0, simStates->mz0, mx1, my1, mz1,
                              demagX, demagY, demagZ, dipoleX, dipoleY, dipoleZ, dmiX, dmiY,
                              dmiZ, t0, simParams->stepsizeHalf, iteration, "1");

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs. This estimates the values of the m-components for the next iteration.
        std::vector<double> mx2(simParams->systemTotalSpins + 2, 0), my2(simParams->systemTotalSpins + 2, 0), mz2(
                simParams->systemTotalSpins + 2, 0);
        RK2StageMultithreaded(mx1, my1, mz1, mx2, my2, mz2,
                              demagX, demagY, demagZ, dipoleX, dipoleY, dipoleZ, dmiX, dmiY,
                              dmiZ, t0, simParams->stepsize, iteration, "2");

        // Everything below heere is part of the method, but not the RK2 stage loops calculations.

        simStates->mx0.clear();
        simStates->my0.clear();
        simStates->mz0.clear();
        mx1.clear();
        my1.clear();
        mz1.clear();

        solverOutputs->SaveDataToFile(mxOutputFile, mx2, iteration);

        // Set final value of current iteration to be starting value of next iteration.
        simStates->mx0 = mx2;
        simStates->my0 = my2;
        simStates->mz0 = mz2;
        mx2.clear();
        my2.clear();
        mz2.clear();

        if ( iteration == simParams->forceStopAtIteration ) {
            std::cout << "Force stop at iteration #" << iteration << std::endl;
            exit(0);
        }

        simParams->totalTime += simParams->stepsize;
    } // End of RK2 FOR loop; all iterations now complete.

    mxOutputFile.close();

    if ( GV.GetEmailWhenCompleted()) { solverOutputs->CreateMetadata(true); }

    if ( simFlags->shouldTrackMagneticMomentNorm ) {
        std::cout << "\nMax norm. value of M is: " << simParams->largestMNorm << std::endl;
    }

    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void SolversImplementation::RK2StageMultithreaded( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                                   const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                                   std::vector<double> &myOut, std::vector<double> &mzOut,
                                                   std::vector<double> &demagX, std::vector<double> &demagY,
                                                   std::vector<double> &demagZ, std::vector<double> &dipoleX,
                                                   std::vector<double> &dipoleY, std::vector<double> &dipoleZ,
                                                   std::vector<double> &dmiX, std::vector<double> &dmiY,
                                                   std::vector<double> &dmiZ, double &currentTime, double &stepsize,
                                                   int &iteration, std::string rkStage ) {

    bool useParallelInvoke = false;  // Mainly for testing purposes at the moment

    if ( simFlags->hasDemagIntense )
        demagField.DemagnetisationFieldIntense(demagX, demagY, demagZ, mxIn, myIn, mzIn);

    std::vector<double> mxInMu, myInMu, mzInMu;
    if ( simFlags->hasDipolar ) {
        mxInMu = mxIn;
        myInMu = myIn;
        mzInMu = mzIn;
        for ( int i = 1; i <= mxIn.size(); i++ ) {
            mxInMu[i] += simParams->PERMITTIVITY_IRON;
            myInMu[i] += simParams->PERMITTIVITY_IRON;
            mzInMu[i] += simParams->PERMITTIVITY_IRON;
        }
    }

    if ( simFlags->hasDMI )
        dmInteraction.calculateOneDimension(mxIn, myIn, mzIn, dmiX, dmiY, dmiZ);


    // Use below line to only use one thread to compare method's results to sequential method
    // tbb::global_control c( tbb::global_control::max_allowed_parallelism, 1 );
    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins),
                      [&]( const tbb::blocked_range<int> r ) {
                          for ( int site = r.begin(); site <= r.end(); site++ ) {
                              // Relative to the current site (site) we have siteLHS and siteRHS
                              int siteLHSLocal = site - 1, siteRHSLocal = site + 1;
                              double dipoleXLocal = 0.0, dipoleYLocal = 0.0, dipoleZLocal = 0.0;
                              double demagXLocal = 0.0, demagYLocal = 0.0, demagZLocal = 0.0;
                              double dmiXLocal = 0.0, dmiYLocal = 0.0, dmiZLocal = 0.0;

                              if ( simFlags->hasDipolar )
                                  dipolarField.DipolarInteractionClassicThreaded(site, mxInMu, myInMu, mzInMu,
                                                                                 dipoleXLocal, dipoleYLocal,
                                                                                 dipoleZLocal);


                              double hxKLocal, hyKLocal, hzKLocal;
                              double mxKLocal, myKLocal, mzKLocal;
                              if ( useParallelInvoke ) {
                                  tbb::parallel_invoke(
                                          // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                                          [&] {
                                              hxKLocal = effectiveField.EffectiveFieldX(site, 0, mxIn[siteLHSLocal],
                                                                                        mxIn[site],
                                                                                        mxIn[siteRHSLocal],
                                                                                        dipoleX[site], demagX[site],
                                                                                        dmiZ[site],
                                                                                        currentTime);
                                          },
                                          [&] {
                                              hyKLocal = effectiveField.EffectiveFieldY(site, 0, myIn[siteLHSLocal],
                                                                                        myIn[site],
                                                                                        myIn[siteRHSLocal],
                                                                                        dipoleY[site], demagY[site],
                                                                                        dmiZ[site]);
                                          },
                                          [&] {
                                              hzKLocal = effectiveField.EffectiveFieldZ(site, 0, mzIn[siteLHSLocal],
                                                                                        mzIn[site],
                                                                                        mzIn[siteRHSLocal],
                                                                                        dipoleZ[site], demagZ[site]);
                                          }
                                  );

                                  tbb::parallel_invoke(
                                          // Calculations for the magnetic moment, coded as symbol 'm', components of the target site
                                          [&] {
                                              mxKLocal = llg.MagneticMomentX(site, mxIn[site], myIn[site], mzIn[site],
                                                                             hxKLocal, hyKLocal,
                                                                             hzKLocal);
                                          },
                                          [&] {
                                              myKLocal = llg.MagneticMomentY(site, mxIn[site], myIn[site], mzIn[site],
                                                                             hxKLocal, hyKLocal,
                                                                             hzKLocal);
                                          },
                                          [&] {
                                              mzKLocal = llg.MagneticMomentZ(site, mxIn[site], myIn[site], mzIn[site],
                                                                             hxKLocal, hyKLocal,
                                                                             hzKLocal);
                                          }
                                  );
                              } else {

                                  // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                                  hxKLocal = effectiveField.EffectiveFieldX(site, 0, mxIn[siteLHSLocal], mxIn[site],
                                                                            mxIn[siteRHSLocal], dipoleXLocal,
                                                                            demagXLocal,
                                                                            dmiZ[site],
                                                                            currentTime);
                                  hyKLocal = effectiveField.EffectiveFieldY(site, 0, myIn[siteLHSLocal], myIn[site],
                                                                            myIn[siteRHSLocal], dipoleYLocal,
                                                                            demagYLocal,
                                                                            dmiZ[site]);
                                  hzKLocal = effectiveField.EffectiveFieldZ(site, 0, mzIn[siteLHSLocal], mzIn[site],
                                                                            mzIn[siteRHSLocal], dipoleZLocal,
                                                                            demagZLocal);

                                  // Calculations for the magnetic moment, coded as symbol 'm', components of the target site
                                  mxKLocal = llg.MagneticMomentX(site, mxIn[site], myIn[site], mzIn[site], hxKLocal,
                                                                 hyKLocal,
                                                                 hzKLocal);
                                  myKLocal = llg.MagneticMomentY(site, mxIn[site], myIn[site], mzIn[site], hxKLocal,
                                                                 hyKLocal,
                                                                 hzKLocal);
                                  mzKLocal = llg.MagneticMomentZ(site, mxIn[site], myIn[site], mzIn[site], hxKLocal,
                                                                 hyKLocal,
                                                                 hzKLocal);
                              }

                              mxOut[site] = simStates->mx0[site] + stepsize * mxKLocal;
                              myOut[site] = simStates->my0[site] + stepsize * myKLocal;
                              mzOut[site] = simStates->mz0[site] + stepsize * mzKLocal;

                              if ( std::isinf(mxOut[site]))
                                  throw std::runtime_error(
                                          "mxOut is INF at site " + std::to_string(site) + " at iteration " +
                                          std::to_string(iteration) +
                                          " in RK2 stage " + rkStage);
                              if ( std::isnan(mxOut[site]))
                                  throw std::runtime_error(
                                          "mxOut is NaN at site " + std::to_string(site) + " at iteration " +
                                          std::to_string(iteration) +
                                          " in RK2 stage " + rkStage);

                              if ( std::isinf(myOut[site]))
                                  throw std::runtime_error(
                                          "myOut is INF at site " + std::to_string(site) + " at iteration " +
                                          std::to_string(iteration) +
                                          " in RK2 stage " + rkStage);
                              if ( std::isnan(myOut[site]))
                                  throw std::runtime_error(
                                          "myOut is NaN at site " + std::to_string(site) + " at iteration " +
                                          std::to_string(iteration) +
                                          " in RK2 stage " + rkStage);

                              if ( std::isinf(mzOut[site]))
                                  throw std::runtime_error(
                                          "mzOut is INF at site " + std::to_string(site) + " at iteration " +
                                          std::to_string(iteration) +
                                          " in RK2 stage " + rkStage);
                              if ( std::isnan(mzOut[site]))
                                  throw std::runtime_error(
                                          "mzOut is NaN at site " + std::to_string(site) + " at iteration " +
                                          std::to_string(iteration) +
                                          " in RK2 stage " + rkStage);
                          }
                      }, tbb::auto_partitioner());

    if ( simFlags->shouldTrackMagneticMomentNorm && rkStage == "2" ) {
        for ( int site = 1; site <= simParams->systemTotalSpins; site++ ) {
            double mIterationNorm = sqrt(pow(mxOut[site], 2) + pow(myOut[site], 2) + pow(mzOut[site], 2));
            if ((simParams->largestMNorm) > (1.0 - mIterationNorm)) {
                simParams->largestMNorm = (1.0 - mIterationNorm);
            }
        }
    }
    // if (simParams->largestMNorm > 1.00005) { throw std::runtime_error("mag. moments are no longer below <= 1.00005"); }

    if ( simFlags->hasDemagFFT )
        std::fill(demagX.begin(), demagX.end(), 0.0);
    std::fill(demagY.begin(), demagY.end(), 0.0);
    std::fill(demagZ.begin(), demagZ.end(), 0.0);
    if ( simFlags->hasDipolar && useParallelInvoke )
        std::fill(dipoleX.begin(), dipoleX.end(), 0.0);
    std::fill(dipoleY.begin(), dipoleY.end(), 0.0);
    std::fill(dipoleZ.begin(), dipoleZ.end(), 0.0);

}

void SolversImplementation::runMethod() {

    std::string methodToUse = GV.GetNumericalMethod();

    if ( methodToUse == "RK2" )
        SolveRK2();
    else if ( methodToUse == "RK2c" )
        SolveRK2Classic();
    else if ( methodToUse == "RK2p" )
        RK2Parallel();
    else
        throw std::runtime_error("Method not recognised");
}