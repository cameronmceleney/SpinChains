//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/SolversImplementation.h"

SolversImplementation::SolversImplementation( std::shared_ptr<SimulationManager> sharedSimManager,
                                              std::shared_ptr<SimulationParameters> sharedSimParams,
                                              std::shared_ptr<SimulationStates> sharedSimStates,
                                              std::shared_ptr<SimulationFlags> sharedSimFlags )

        : SolversSuperClass(std::move(sharedSimManager), std::move(sharedSimParams), std::move(sharedSimStates), std::move(sharedSimFlags)),
          demagField(simParams.get(), simStates.get(), simFlags.get()),
          dmInteraction(simParams.get(), simStates.get(), simFlags.get()),
          effectiveField(simParams.get(), simStates.get(), simFlags.get()),
          dipolarField(simParams.get(), simStates.get(), simFlags.get()),
          llg(simParams.get(), simStates.get(), simFlags.get()),
          stt(simParams.get(), simStates.get(), simFlags.get()),
          exchangeField(simParams.get(), simStates.get(), simFlags.get()),
          biasField(simParams.get(), simStates.get(), simFlags.get()) {}

void SolversImplementation::_testShockwaveConditions( double iteration ) {

    if ( simFlags->shouldDriveCeaseEarly ) {
        // and (simParams->isShockwaveOn and simParams->isShockwaveAtMax)) {
        if ( simFlags->isShockwaveOn and not simFlags->isShockwaveAtMax )
            std::cout << "Shock not at maximum when cut-off" << std::endl;

        if ( iteration >= simParams->iterationEnd * simParams->risingTimeEndAtIteration ) {
            // Shockwave ends once simulation is a certain % complete
            simFlags->hasRisingTime = false;
            simFlags->isShockwaveOn = false;
            simParams->oscillatingZeemanStrength = 0;
        }

        return;
    }

    // If method is triggered, then the applied biasFieldDriving is increased by the scale factor risingTimeScalingFactor
    if (simFlags->hasRisingTime and not simFlags->isShockwaveOn ) {
        if ( iteration >= simParams->iterationEnd * simParams->risingTimeStartAtIteration ) {
            // Shockwave begins once simulation is a certain % complete
            simFlags->isShockwaveOn = true;
            simParams->oscillatingZeemanStrength = simParams->risingTimeInitialMagnitude;
        }

        return;
    }

    if ( simFlags->isShockwaveOn and not simFlags->isShockwaveAtMax ) {
        simParams->oscillatingZeemanStrength += simParams->shockwaveStepsize;

        if ( simParams->oscillatingZeemanStrength >= simParams->risingTimeMaximum ) {
            simParams->oscillatingZeemanStrength = simParams->risingTimeMaximum;
            simFlags->isShockwaveAtMax = true;
        }
        return;

    }

}

void SolversImplementation::SolveRK2Classic() {
    std::shared_ptr<SolversDataHandling> childNMData = std::make_shared<SolversDataHandling>(simManager, simParams, simStates,
                                                                                             simFlags);
    // Only uses a single spin chain to solve the RK2 midpoint method.
    CommonStructures::Timer methodTimer;

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

    methodTimer.setName("RK2 Sequential (Classic)");
    methodTimer.start();

    std::vector<double> demagX(simParams->systemTotalSpins + 2, 0.0), demagY(simParams->systemTotalSpins + 2, 0.0), demagZ(
            simParams->systemTotalSpins + 2, 0.0);
    //std::vector<double> dipoleX(simParams->systemTotalSpins + 2, 0.0), dipoleY(simParams->systemTotalSpins + 2, 0.0), dipoleZ(simParams->systemTotalSpins + 2, 0.0);

    for ( int iteration = simParams->iterationStart; iteration <= simParams->iterationEnd; iteration++ ) {

        if ( simParams->iterationEnd >= 100 && iteration % (simParams->iterationEnd / 100) == 0 )
            bar.update(); // Doesn't work for fewer than 100 iterations

        _testShockwaveConditions(iteration);

        double t0 = simParams->totalTime;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = simParams->mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(simParams->systemTotalSpins + 2, 0), my1(simParams->systemTotalSpins + 2, 0), mz1(simParams->systemTotalSpins + 2, 0);
        if ( simFlags->hasDemagFactors ) {
            demagField.calculateOneDimension(simStates->mx0, simStates->my0, simStates->mz0, demagX, demagY, demagZ,
                                             DemagnetisationFields::DemagMethod::Demag,
                                             CommonStructures::Parallelisations::Sequential);
        } else if ( simFlags->hasDemagFFT ) {
            std::string rkStageName = "2-1";
            demagField.calculateOneDimension(simStates->mx0, simStates->my0, simStates->mz0, demagX, demagY, demagZ,
                                             DemagnetisationFields::DemagMethod::DipoleGreenFunctionReal,
                                             CommonStructures::Parallelisations::Sequential);        }

        // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)
        // RK2 Stage 1. Takes initial conditions as inputs.
        for ( int site = 1; site <= simParams->systemTotalSpins; site++ ) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

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
                CommonStructures::Vector2D mxTermsForDMI = {simStates->mx0[spinLHS], simStates->mx0[site]};
                CommonStructures::Vector2D myTermsForDMI = {simStates->my0[spinLHS], simStates->my0[site]};
                CommonStructures::Vector2D mzTermsForDMI = {simStates->mz0[spinLHS], simStates->mz0[site]};
                CommonStructures::Vector3D dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                                myTermsForDMI,
                                                                                mzTermsForDMI);

                dmiZ = dmiTerms.z();
            }

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK0 = effectiveField.EffectiveFieldXClassic(site, 0, simStates->mx0[spinLHS], simStates->mx0[site],
                                                                simStates->mx0[spinRHS], dipoleX, demagX[site],
                                                                dmiZ, t0);
            double hyK0 = effectiveField.EffectiveFieldYClassic(site, 0, simStates->my0[spinLHS], simStates->my0[site],
                                                                simStates->my0[spinRHS], dipoleY, demagY[site],
                                                                dmiZ);
            double hzK0 = effectiveField.EffectiveFieldZClassic(site, 0, simStates->mz0[spinLHS], simStates->mz0[site],
                                                                simStates->mz0[spinRHS], dipoleZ, demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double localGilbertFactor = llg._checkIfDampingMapExists(site);
            double mxK1 = llg.MagneticMomentX(site, simStates->mx0[site], simStates->my0[site], simStates->mz0[site],
                                              hxK0, hyK0, hzK0, localGilbertFactor);
            double myK1 = llg.MagneticMomentY(site, simStates->mx0[site], simStates->my0[site], simStates->mz0[site],
                                              hxK0, hyK0, hzK0, localGilbertFactor);
            double mzK1 = llg.MagneticMomentZ(site, simStates->mx0[site], simStates->my0[site], simStates->mz0[site],
                                              hxK0, hyK0, hzK0, localGilbertFactor);

            // Find (m0 + k1/2) for each site, which is used in the next stage.
            mx1[site] = simStates->mx0[site] + simParams->stepsizeHalf * mxK1;
            my1[site] = simStates->my0[site] + simParams->stepsizeHalf * myK1;
            mz1[site] = simStates->mz0[site] + simParams->stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        // EASY FIND
        std::vector<double> mx2(simParams->systemTotalSpins + 2, 0), my2(simParams->systemTotalSpins + 2, 0), mz2(simParams->systemTotalSpins + 2, 0);
        std::fill(demagX.begin(), demagX.end(), 0.0);
        std::fill(demagY.begin(), demagY.end(), 0.0);
        std::fill(demagZ.begin(), demagZ.end(), 0.0);
        //std::fill(dipoleX.begin(), dipoleX.end(), 0.0); std::fill(dipoleY.begin(), dipoleY.end(), 0.0); std::fill(dipoleZ.begin(), dipoleZ.end(), 0.0);

        //if (simFlags->hasDipolar)
        //    dipolarField.DipolarInteraction1D(mx1, dipoleX);
        if ( simFlags->hasDemagFactors ) {
            demagField.calculateOneDimension(simStates->mx0, simStates->my0, simStates->mz0, demagX, demagY, demagZ,
                                             DemagnetisationFields::DemagMethod::Demag,
                                             CommonStructures::Parallelisations::Sequential);
        } else if ( simFlags->hasDemagFFT ) {
            std::string rkStageName = "2-2";
            demagField.calculateOneDimension(simStates->mx0, simStates->my0, simStates->mz0, demagX, demagY, demagZ,
                                             DemagnetisationFields::DemagMethod::DipoleGreenFunctionReal,
                                             CommonStructures::Parallelisations::Sequential);
            // if (iteration > 0) {std::cout << "Stage 2" << std::endl; PrintVector(demagZ, false);}
        }

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
        for ( int site = 1; site <= simParams->systemTotalSpins; site++ ) {

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
                CommonStructures::Vector2D mxTermsForDMI = {simStates->mx0[spinLHS], simStates->mx0[site]};
                CommonStructures::Vector2D myTermsForDMI = {simStates->my0[spinLHS], simStates->my0[site]};
                CommonStructures::Vector2D mzTermsForDMI = {simStates->mz0[spinLHS], simStates->mz0[site]};
                CommonStructures::Vector3D dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                                myTermsForDMI,
                                                                                mzTermsForDMI);

                dmiZ = dmiTerms[2];
            }

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK1 = effectiveField.EffectiveFieldXClassic(site, 0, mx1[spinLHS], mx1[site], mx1[spinRHS], dipoleX,
                                                                demagX[site], dmiZ, t0);
            double hyK1 = effectiveField.EffectiveFieldYClassic(site, 0, my1[spinLHS], my1[site], my1[spinRHS], dipoleY,
                                                                demagY[site], dmiZ);
            double hzK1 = effectiveField.EffectiveFieldZClassic(site, 0, mz1[spinLHS], mz1[site], mz1[spinRHS], dipoleZ,
                                                                demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double localGilbertFactor = llg._checkIfDampingMapExists(site);
            double mxK2 = llg.MagneticMomentX(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1, localGilbertFactor);
            double myK2 = llg.MagneticMomentY(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1, localGilbertFactor);
            double mzK2 = llg.MagneticMomentZ(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1, localGilbertFactor);

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
    methodTimer.stop();

    if ( GV.GetEmailWhenCompleted()) {
        childNMData->CreateMetadata(true);
    }

    if ( simFlags->shouldTrackMagneticMomentNorm )
        std::cout << "\nMax norm. value of M is: " << simParams->largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
    methodTimer.print();
}

void SolversImplementation::SolveRK2() {
    std::shared_ptr<SolversDataHandling> childNMData = std::make_shared<SolversDataHandling>(simManager, simParams, simStates,
                                                                                             simFlags);
    // Uses multiple layers to solve the RK2 midpoint method. See the documentation for more details.
    CommonStructures::Timer rk2Timer;

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

    rk2Timer.setName("RK2 Sequential (Layers)");
    rk2Timer.start();

    std::vector<double> demagX(simParams->systemTotalSpins, 0.0), demagY(simParams->systemTotalSpins, 0.0), demagZ(
            simParams->systemTotalSpins + 2, 0.0);

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
                    CommonStructures::Vector2D mxTermsForDMI = {simStates->m0Nest[layer][spinLHS][0],
                                                           simStates->m0Nest[layer][site][0]};
                    CommonStructures::Vector2D myTermsForDMI = {simStates->m0Nest[layer][spinLHS][1],
                                                           simStates->m0Nest[layer][site][1]};
                    CommonStructures::Vector2D mzTermsForDMI = {simStates->m0Nest[layer][spinLHS][2],
                                                           simStates->m0Nest[layer][site][2]};
                    CommonStructures::Vector3D dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                                    myTermsForDMI,
                                                                                    mzTermsForDMI);

                    dmiZ = dmiTerms[2];
                }

                // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                double hxK0 = effectiveField.EffectiveFieldXClassic(site, layer, mxLHS, mxMID, mxRHS, dipoleX,
                                                                    demagX[site], dmiZ, t0);
                double hyK0 = effectiveField.EffectiveFieldYClassic(site, layer, myLHS, myMID, myRHS, dipoleY,
                                                                    demagY[site], dmiZ);
                double hzK0 = effectiveField.EffectiveFieldZClassic(site, layer, mzLHS, mzMID, mzRHS, dipoleZ,
                                                                    demagZ[site]);

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
                    CommonStructures::Vector2D mxTermsForDMI = {simStates->m0Nest[layer][spinLHS][0],
                                                           simStates->m0Nest[layer][site][0]};
                    CommonStructures::Vector2D myTermsForDMI = {simStates->m0Nest[layer][spinLHS][1],
                                                           simStates->m0Nest[layer][site][1]};
                    CommonStructures::Vector2D mzTermsForDMI = {simStates->m0Nest[layer][spinLHS][2],
                                                           simStates->m0Nest[layer][site][2]};
                    CommonStructures::Vector3D dmiTerms = dmInteraction.calculateClassic(site, mxTermsForDMI,
                                                                                    myTermsForDMI,
                                                                                    mzTermsForDMI);

                    dmiZ = dmiTerms[2];
                }
                // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                double hxK1 = effectiveField.EffectiveFieldXClassic(site, layer, mxLHS, mxMID, mxRHS, dipoleX,
                                                                    demagX[site], dmiZ, t0);
                double hyK1 = effectiveField.EffectiveFieldYClassic(site, layer, myLHS, myMID, myRHS, dipoleY,
                                                                    demagY[site], dmiZ);
                double hzK1 = effectiveField.EffectiveFieldZClassic(site, layer, mzLHS, mzMID, mzRHS, dipoleZ,
                                                                    demagZ[site]);

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

    rk2Timer.stop();

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
    rk2Timer.print();
}

void SolversImplementation::RK2Parallel() {
    CommonStructures::Timer parallelTimer;
    // Only works for a 1D spin chain
    std::shared_ptr<SolversDataHandling> solverOutputs = std::make_shared<SolversDataHandling>(simManager, simParams, simStates,
                                                                                               simFlags);

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxOutputFile(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");

    if ( simManager->massProduce || simFlags->resetSimState )
        _clearDataStorageVariables();

    if ( simFlags->isFerromagnetic ) {
        if ( !simManager->massProduce ) { solverOutputs->InformUserOfCodeType("RK2 Midpoint (Parallel)(FM)"); }
        solverOutputs->CreateFileHeader(mxOutputFile, "RK2 Midpoint (Parallel)(FM)");
    } else {
        if ( !simManager->massProduce ) { solverOutputs->InformUserOfCodeType("RK2 Midpoint (AFM)"); }
        solverOutputs->CreateFileHeader(mxOutputFile, "RK2 Midpoint (AFM)");
    }

    if ( GV.GetEmailWhenCompleted()) { solverOutputs->CreateMetadata(); }

    if ( !simManager->massProduce ) {
        simProgressBar.reset();
        simProgressBar.set_niter(progressBarSubdivisions);
        simProgressBar.show_bar();
    }

    parallelTimer.setName("RK2 Midpoint (Parallel)");
    parallelTimer.start();

    // TODO. As this is all multi-threaded, having such larger vectors is excessive. Change methods to work by element (so that onl mx/my/mz are large in memory)
    // Only create vectors once to reuse memory. Only assign memory if flags are true. Faster to declare loop-wide vectors than define them in each loop iteration.
    bool useCompact = true;
    if ( useCompact )
        _resizeClassContainersTest();
    else
        _resizeClassContainers();

    demagField.initialise(false);

    for ( int iteration = simParams->iterationStart; iteration <= simParams->iterationEnd; iteration++ ) {

        if ( !simManager->massProduce && simParams->iterationEnd >= progressBarSubdivisions &&
             iteration % (simParams->iterationEnd / progressBarSubdivisions) == 0 )
            simProgressBar.update(); // Need to match definition of bar

        _testShockwaveConditions(iteration);

        double t0 = simParams->totalTime, t0Half = simParams->totalTime + simParams->stepsizeHalf;

        if ( useCompact ) {
            // RK2 Stage 1. Takes initial conditions as inputs. The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = simParams->mx0 + (h * k1 / 2) etc
            // Possible mistake here; should it not be 't0+h' and not 't0' for RK-S1?
            RK2StageMultithreadedCompact(simStates->mx0, simStates->my0, simStates->mz0, mx1p, my1p, mz1p,
                                         t0, simParams->stepsizeHalf, iteration, "1");

            // RK2 Stage 2. Takes (m0 + k1/2) as inputs. This estimates the values of the m-components for the next iteration.
            RK2StageMultithreadedCompact(mx1p, my1p, mz1p, mx2p, my2p, mz2p,
                                         t0Half, simParams->stepsize, iteration, "2");
        } else {
            // RK2 Stage 1. Takes initial conditions as inputs. The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = simParams->mx0 + (h * k1 / 2) etc
            RK2StageMultithreaded(simStates->mx0, simStates->my0, simStates->mz0, mx1p, my1p, mz1p,
                                  simParams->totalTime, simParams->stepsizeHalf, iteration, "1");
            // RK2 Stage 2. Takes (m0 + k1/2) as inputs. This estimates the values of the m-components for the next iteration.
            RK2StageMultithreaded(mx1p, my1p, mz1p, mx2p, my2p, mz2p,
                                  simParams->totalTime, simParams->stepsize, iteration, "2");
        }

        // Everything below here is part of the method, but not the RK2 stage loops calculations.

        // Do not clear mx0/mx1/mx2 (etc) if they are never resized and only to be refilled.

        solverOutputs->SaveDataToFile(mxOutputFile, mx2p, iteration);

        // Set final value of current iteration to be starting value of next iteration.
        simStates->mx0 = mx2p;
        simStates->my0 = my2p;
        simStates->mz0 = mz2p;

        if ( iteration == simParams->forceStopAtIteration ) {
            std::cout << "Force stop at iteration #" << iteration << std::endl;
            exit(0);
        }

        simParams->totalTime += simParams->stepsize;
    } // End of RK2 FOR loop; all iterations now complete.

    mxOutputFile.close();

    parallelTimer.stop();

    /* Should turn this code into a debug method
    for (int i = 1; i <= simParams->systemTotalSpins; i++) {
        auto it = simStates->dRGradientMap.find(i);
        if (it != simStates->dRGradientMap.end()) {
            // We found this site in the map and in the driving region
            std::cout << "Site " << i << " has a gradient of " << it->second << std::endl;
        }
    }
     */

    if ( GV.GetEmailWhenCompleted()) { solverOutputs->CreateMetadata(true); }

    if ( !simManager->massProduce ) {
        if ( simFlags->shouldTrackMagneticMomentNorm ) {
            std::cout << "\nMax norm. value of M is: " << simParams->largestMNorm << std::endl;
        }

        std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;

        parallelTimer.print();
    }
}

void SolversImplementation::RK4Parallel() {
    CommonStructures::Timer parallelTimer;
    // Only works for a 1D spin chain
    std::shared_ptr<SolversDataHandling> solverOutputs = std::make_shared<SolversDataHandling>(simManager, simParams, simStates,
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

    if ( !simManager->massProduce ) {
        simProgressBar.reset();
        simProgressBar.set_niter(progressBarSubdivisions);
        simProgressBar.show_bar();
    };

    parallelTimer.setName("RK4 Midpoint (Parallel)");
    parallelTimer.start();

    // TODO. As this is all multi-threaded, having such larger vectors is excessive. Change methods to work by element (so that onl mx/my/mz are large in memory)
    // Only create vectors once to reuse memory. Only assign memory if flags are true. Faster to declare loop-wide vectors than define them in each loop iteration.
    _resizeClassContainersRK4();

    for ( int iteration = simParams->iterationStart; iteration <= simParams->iterationEnd; iteration++ ) {

        if ( !simManager->massProduce && simParams->iterationEnd >= progressBarSubdivisions && iteration % (simParams->iterationEnd / progressBarSubdivisions) == 0 )
            simProgressBar.update(); // Doesn't work for fewer than 100 iterations

        _testShockwaveConditions(iteration);

        double stepsizeZero = 0;

        double t0 = simParams->totalTime, t0Half = simParams->totalTime + simParams->stepsizeHalf, t0Full =
                simParams->totalTime + simParams->stepsize;

        // RK4 Stage 1. Takes initial conditions as inputs. The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = simParams->mx0 + (h * k1 / 2) etc
        // Possible mistake here; should it not be 't0+h' and not 't0' for RK-S1?
        RK4StageMultithreadedCompact(simStates->mx0, simStates->my0, simStates->mz0, mx1p, my1p, mz1p,
                                     t0, simParams->stepsize, iteration, "1");

        // RK2 Stage 4. Takes (m0 + k2/2) as inputs. This estimates the values of the m-components for the next iteration.
        RK4StageMultithreadedCompact(mx1p, my1p, mz1p, mx2p, my2p, mz2p,
                                     t0Half, simParams->stepsize, iteration, "2");

        // RK2 Stage 3. Takes (m0 + k3/2) as inputs
        RK4StageMultithreadedCompact(mx2p, my2p, mz2p, mx3p, my3p, mz3p,
                                     t0Half, simParams->stepsize, iteration, "3");

        // RK2 Stage 4. Takes (m0 + k3) as inputs
        RK4StageMultithreadedCompact(mx3p, my3p, mz3p, mx4p, my4p, mz4p,
                                     t0Full, simParams->stepsize, iteration, "4");

        // Everything below here is part of the method, but not the RK4 stage loops calculations.

        // Do not clear mx0/mx1/mx2 (etc) if they are never resized and only to be refilled.

        solverOutputs->SaveDataToFile(mxOutputFile, mx2p, iteration);

        // Set final value of current iteration to be starting value of next iteration.
        simStates->mx0 = mx2p;
        simStates->my0 = my2p;
        simStates->mz0 = mz2p;

        if ( iteration == simParams->forceStopAtIteration ) {
            std::cout << "Force stop at iteration #" << iteration << std::endl;
            exit(0);
        }

        simParams->totalTime += simParams->stepsize;
    } // End of RK2 FOR loop; all iterations now complete.

    mxOutputFile.close();

    parallelTimer.stop();

    if ( GV.GetEmailWhenCompleted()) { solverOutputs->CreateMetadata(true); }

    if ( simFlags->shouldTrackMagneticMomentNorm ) {
        std::cout << "\nMax norm. value of M is: " << simParams->largestMNorm << std::endl;
    }

    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
    parallelTimer.print();
}

void SolversImplementation::RK2StageMultithreaded( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                                   const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                                   std::vector<double> &myOut, std::vector<double> &mzOut,
                                                   double &currentTime, double &stepsize, int &iteration,
                                                   std::string rkStage ) {
    CommonStructures::Timer dipolarTimer;
    bool useParallel = true;

    if ( simFlags->hasDemagFactors )
        demagField.calculateOneDimension(simStates->mx0, simStates->my0, simStates->mz0, demagXp, demagYp, demagZp,
                                 DemagnetisationFields::DemagMethod::Demag,
                                 CommonStructures::Parallelisations::Multithreaded);
    std::vector<double> mxInMu, myInMu, mzInMu;
    if ( simFlags->hasDipolar ) {
        // Required step to convert the magnetic moment components to the magnetic field components
        //mxInMu = mxIn;
        //myInMu = myIn;
        //mzInMu = mzIn;
        //for ( int i = 1; i <= mxIn.size(); i++ ) {
        //    mxInMu[i] *= simParams.PERMEABILITY_IRON;//simParams->satMag * simParams->systemTotalSpins;
        //    myInMu[i] *= simParams.PERMEABILITY_IRON;//simParams->satMag * simParams->systemTotalSpins;
        //    mzInMu[i] *= simParams.PERMEABILITY_IRON;//simParams->satMag * simParams->systemTotalSpins;
        //}
        /// Option 1
        dipolarField.DipolarInteractionClassicThreaded(mxIn, myIn, mzIn, dipoleXp, dipoleYp, dipoleZp);
    }

    if ( simFlags->hasDMI )
        dmInteraction.calculateOneDimension(mxIn, myIn, mzIn, dmiXp, dmiYp, dmiZp, useParallel);

    if ( simFlags->hasSTT )
        // placeholder for example. Find STT for each site at all sites before main loop
        stt.calculateOneDimension(mxIn, myIn, mzIn, sttXp, sttYp, sttZp);

    dipolarTimer.setName("Dipolar");

    // Use below line to only use one thread to compare method's results to sequential method
    // tbb::global_control c( tbb::global_control::max_allowed_parallelism, 1 );
    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins), [&]( const tbb::blocked_range<int> r ) {
        for ( int site = r.begin(); site <= r.end(); site++ ) {
            /*
            * All declarations of overwritten variables are placed here as 'local' to ensure
            * that each thread has its own definition; else variables are not threadsafe
            */

            // Relative to the current site (site) we have siteLHS and siteRHS
            int siteLHSLocal = site - 1, siteRHSLocal = site + 1;
            DipoleTerms dipoleLocal;
            DemagTerms demagLocal;
            DMITerms dmiLocal;
            STTTerms sttLocal;

            // All need to be defined as default cases in case their flags aren't called to overwrite
            if ( simFlags->hasDipolar ) {
                dipoleLocal.x = dipoleXp[site];
                dipoleLocal.y = dipoleYp[site];
                dipoleLocal.z = dipoleZp[site];
            }
            if ( simFlags->hasDemagFactors ) {
                demagLocal.x = demagXp[site];
                demagLocal.y = demagYp[site];
                demagLocal.z = demagZp[site];
            }
            if ( simFlags->hasDMI ) { dmiLocal.z = dmiZp[site]; }
            if ( simFlags->hasSTT ) {
                sttLocal.x = sttXp[site];
                sttLocal.y = sttYp[site];
                sttLocal.z = sttZp[site];
            }

            // Will always be initialised so only need to declare here
            HkTerms hkLocal;
            MkTerms mkLocal;

            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            hkLocal.x = effectiveField.EffectiveFieldXClassic(site, 0, mxIn[siteLHSLocal], mxIn[site],
                                                              mxIn[siteRHSLocal], dipoleLocal.x, demagLocal.x,
                                                              dmiLocal.z, currentTime);
            hkLocal.y = effectiveField.EffectiveFieldYClassic(site, 0, myIn[siteLHSLocal], myIn[site],
                                                              myIn[siteRHSLocal], dipoleLocal.y, demagLocal.y,
                                                              dmiLocal.z);
            hkLocal.z = effectiveField.EffectiveFieldZClassic(site, 0, mzIn[siteLHSLocal], mzIn[site],
                                                              mzIn[siteRHSLocal], dipoleLocal.z, demagLocal.z);

            // Calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double localGilbertFactor = llg._checkIfDampingMapExists(site);
            mkLocal.x = llg.MagneticMomentX(site, mxIn[site], myIn[site], mzIn[site],
                                            hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);
            mkLocal.y = llg.MagneticMomentY(site, mxIn[site], myIn[site], mzIn[site],
                                            hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);
            mkLocal.z = llg.MagneticMomentZ(site, mxIn[site], myIn[site], mzIn[site],
                                            hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);


            mxOut[site] = simStates->mx0[site] + stepsize * mkLocal.x;
            myOut[site] = simStates->my0[site] + stepsize * mkLocal.y;
            mzOut[site] = simStates->mz0[site] + stepsize * mkLocal.z;

            _testOutputValues(mxOut[site], myOut[site], mzOut[site], site, iteration, rkStage);
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

    // No need to fill demagX/dipoleX (etc) here they are constantly overwritten; filling to zeroes is just a waste unless debugging

}

void SolversImplementation::RK2StageMultithreadedTest( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                                       const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                                       std::vector<double> &myOut, std::vector<double> &mzOut,
                                                       const double &currentTime, const double &stepsize,
                                                       const int &iteration, std::string rkStage ) {
    CommonStructures::Timer dipolarTimer;
    bool useParallel = true;
    int layer = 0;

    // Use these to ensure that there's no issues with sharing a member attribute across
    // methods (this is less memory efficient, but for debugging purposes)
    //std::vector<std::atomic<double>> effectiveFieldXAtomic(simParams->systemTotalSpins + 2);
    //std::vector<std::atomic<double>> effectiveFieldYAtomic(simParams->systemTotalSpins + 2);
    //std::vector<std::atomic<double>> effectiveFieldZAtomic(simParams->systemTotalSpins + 2);

    //tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins), [&]( const tbb::blocked_range<int> tbbRange ) {
    //    for ( int i = tbbRange.begin(); i <= tbbRange.end(); i++ ) {
//
    //        CommonStructures::Vector3D driveTemp = effectiveField.EffectiveFieldsCombinedTestDriveOnly(i, 0, mxIn, myIn,
    //                                                                                              mzIn, currentTime);
    //        effectiveFieldXLocal[i] = driveTemp[0];
    //        effectiveFieldYLocal[i] = driveTemp[1];
    //        effectiveFieldZLocal[i] = driveTemp[2];[image] ([pdf](zotero://open-pdf/library/items/6SHVBE95?page=3&annotation=95Z2LQSE))
    //([Moon et al., 2013, p. 3](zotero://select/library/items/A7DNR44Z))
    //    }
    //}, tbb::auto_partitioner());
    //exchangeField.calculateOneDimension(mxIn, myIn, mzIn, effectiveFieldXLocal, effectiveFieldYLocal, effectiveFieldZLocal, useParallel);



    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins),
                      [&]( const tbb::blocked_range<int> tbbRange ) {
                          for ( int i = tbbRange.begin(); i <= tbbRange.end(); i++ ) {

                              CommonStructures::Vector3D driveTemp = effectiveField.EffectiveFieldsCombinedTestDriveOnly(i,
                                                                                                                    0,
                                                                                                                    mxIn,
                                                                                                                    myIn,
                                                                                                                    mzIn,
                                                                                                                    currentTime);
                              CommonStructures::Vector3D exchangeTemp = effectiveField.EffectiveFieldsCombinedTestExOnly(i,
                                                                                                                    0,
                                                                                                                    mxIn,
                                                                                                                    myIn,
                                                                                                                    mzIn);
                              effectiveFieldX[i] = driveTemp[0] + exchangeTemp[0];
                              effectiveFieldY[i] = driveTemp[1] + exchangeTemp[1];
                              effectiveFieldZ[i] = driveTemp[2] + exchangeTemp[2];
                          }
                      }, tbb::auto_partitioner());


    /*
    std::vector<std::atomic<double>> effectiveFieldXLocal(simParams->systemTotalSpins);
    std::vector<std::atomic<double>> effectiveFieldYLocal(simParams->systemTotalSpins);
    std::vector<std::atomic<double>> effectiveFieldZLocal(simParams->systemTotalSpins);

    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            // All computations where calculations are site-by-site should be here
            CommonStructures::Vector3D exchangeTemp = effectiveField.EffectiveFieldsCombinedTestExOnly(i, 0, mxIn, myIn, mzIn);

            effectiveFieldXLocal[i].fetch_add(exchangeTemp[0]);
            effectiveFieldYLocal[i].fetch_add(exchangeTemp[1]);
            effectiveFieldZLocal[i].fetch_add(exchangeTemp[2]);
        }
    }, tbb::auto_partitioner());

    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            // All computations where calculations are site-by-site should be here
            CommonStructures::Vector3D driveTemp = effectiveField.EffectiveFieldsCombinedTestDriveOnly(i, 0, mxIn, myIn, mzIn, currentTime);

            effectiveFieldXLocal[i].fetch_add(driveTemp[0]);
            effectiveFieldYLocal[i].fetch_add(driveTemp[1]);
            effectiveFieldZLocal[i].fetch_add(driveTemp[2]);
        }
    }, tbb::auto_partitioner());

    _transferDataThenReleaseAtomicVector(effectiveFieldXLocal, effectiveFieldX);
    _transferDataThenReleaseAtomicVector(effectiveFieldYLocal, effectiveFieldY);
    _transferDataThenReleaseAtomicVector(effectiveFieldZLocal, effectiveFieldZ);
     */

    /*
    tbb::parallel_reduce(tbb::blocked_range<int>(1, simParams->systemTotalSpins), 0,
        [&](const tbb::blocked_range<int>& range, int dummy) -> int {
            for (int i = range.begin(); i < range.end(); ++i) {
                CommonStructures::Vector3D exchangeTemp = effectiveField.EffectiveFieldsCombinedTestExOnly(i, 0, mxIn, myIn, mzIn);
                CommonStructures::Vector3D driveTemp = effectiveField.EffectiveFieldsCombinedTestDriveOnly(i, 0, mxIn, myIn, mzIn, currentTime);

                effectiveFieldX[i] += exchangeTemp[0] + driveTemp[0];
                effectiveFieldY[i] += exchangeTemp[1] + driveTemp[1];
                effectiveFieldZ[i] += exchangeTemp[2] + driveTemp[2];
            }
            return 0; // The dummy variable remains unchanged
        },
        [](int x, int y) -> int {
            return 0; // The dummy variable remains unchanged
        },
        tbb::auto_partitioner()
    );
     */
    /*
    exchangeField.calculateOneDimension(mxIn, myIn, mzIn, effectiveFieldXAtomic, effectiveFieldYAtomic, effectiveFieldZAtomic, useParallel);
    biasField.calculateOneDimension(layer, currentTime, mzIn, effectiveFieldXAtomic, effectiveFieldYAtomic,
                                    effectiveFieldZAtomic, useParallel);

    std::vector<double> effectiveFieldX(simParams->systemTotalSpins+2);
    std::vector<double> effectiveFieldY(simParams->systemTotalSpins+2);
    std::vector<double> effectiveFieldZ(simParams->systemTotalSpins+2);
    _transferDataThenReleaseAtomicVector(effectiveFieldXAtomic, effectiveFieldX, true);
    _transferDataThenReleaseAtomicVector(effectiveFieldYAtomic, effectiveFieldY, true);
    _transferDataThenReleaseAtomicVector(effectiveFieldZAtomic, effectiveFieldZ, true);
     */

    dipolarTimer.setName("Dipolar");

    // Testing. Probably should use single thread, as the overhead here likely won't be worthwhile
    //tbb::global_control c( tbb::global_control::max_allowed_parallelism, 1 );
    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins),
                      [&]( const tbb::blocked_range<int> tbbRange ) {
                          for ( int site = tbbRange.begin(); site <= tbbRange.end(); site++ ) {
                              /*
                               * All declarations of overwritten variables are placed here as 'local' to ensure
                               * that each thread has its own definition; else variables are not threadsafe
                              */

                              // Will always be initialised so only need to declare here
                              //HkTerms hkLocal;
                              //MkTerms mkLocal;

                              double hkLocalX = effectiveFieldX[site];
                              double hkLocalY = effectiveFieldY[site];
                              double hkLocalZ = effectiveFieldZ[site];

                              // Calculations for the magnetic moment, coded as symbol 'm', components of the target site
                              double localGilbertFactor = llg._checkIfDampingMapExists(site);
                              double mkLocalX = llg.MagneticMomentX(site, mxIn[site], myIn[site], mzIn[site], hkLocalX,
                                                                    hkLocalY, hkLocalZ, localGilbertFactor);
                              double mkLocalY = llg.MagneticMomentY(site, mxIn[site], myIn[site], mzIn[site], hkLocalX,
                                                                    hkLocalY, hkLocalZ, localGilbertFactor);
                              double mkLocalZ = llg.MagneticMomentZ(site, mxIn[site], myIn[site], mzIn[site], hkLocalX,
                                                                    hkLocalY, hkLocalZ, localGilbertFactor);

                              mxOut[site] = simStates->mx0[site] + stepsize * mkLocalX;
                              myOut[site] = simStates->my0[site] + stepsize * mkLocalY;
                              mzOut[site] = simStates->mz0[site] + stepsize * mkLocalZ;

                              _testOutputValues(mxOut[site], myOut[site], mzOut[site], site, iteration, rkStage);
                          }
                      }, tbb::auto_partitioner());

    if ( !simFlags->shouldTrackMagneticMomentNorm && rkStage == "2" ) {
        for ( int site = 1; site <= simParams->systemTotalSpins; site++ ) {
            double mIterationNorm = sqrt(pow(mxOut[site], 2) + pow(myOut[site], 2) + pow(mzOut[site], 2));
            if ((simParams->largestMNorm) > (1.0 - mIterationNorm)) {
                simParams->largestMNorm = (1.0 - mIterationNorm);
            }
        }
    }

    // if (simParams->largestMNorm > 1.00005) { throw std::runtime_error("mag. moments are no longer below <= 1.00005"); }

    // No need to fill demagX/dipoleX (etc) here they are constantly overwritten; filling to zeroes is just a waste unless debugging

}

void
SolversImplementation::RK2StageMultithreadedCompact( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                                     const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                                     std::vector<double> &myOut, std::vector<double> &mzOut,
                                                     const double &currentTime, const double &stepsize,
                                                     const int &iteration, std::string rkStage ) {
    CommonStructures::Timer dipolarTimer;
    bool useParallel = true;
    int layer = 0;

    // Use these to ensure that there's no issues with sharing a member attribute across
    // methods (this is less memory efficient, but for debugging purposes)
    // std::vector<std::atomic<double>> effectiveFieldXAtomic(simParams->systemTotalSpins + 2);
    // std::vector<std::atomic<double>> effectiveFieldYAtomic(simParams->systemTotalSpins + 2);
    // std::vector<std::atomic<double>> effectiveFieldZAtomic(simParams->systemTotalSpins + 2);
    std::vector<double> dipoleFieldX, dipoleFieldY, dipoleFieldZ;
    if (simFlags->hasDemagFFT) {
        std::tie(dipoleFieldX, dipoleFieldY, dipoleFieldZ) = demagField.calculateDipoleFieldGreenReal(1, 1, mxIn, myIn,mzIn);
    }
    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins + 1),
            [&]( const tbb::blocked_range<int> tbbRange ) {
                for ( int site = tbbRange.begin(); site < tbbRange.end(); site++ ) {
                    auto exchangeTemp = exchangeField.calculateExchangeField(1, site, mxIn, myIn, mzIn);
                    auto biasTemp = biasField.calculateBiasField(1, site, layer, currentTime, mzIn);

                    CommonStructures::Vector3D dipolarTemp{0.0, 0.0, 0.0};
                    if (simFlags->hasDemagFactors) {
                        // placeholder as demagnetisation-field code hasn't been adapted for direct invocation
                        dipolarTemp = demagField.calculateDemagField(1, site, mxIn, myIn, mzIn);
                    } else if (simFlags->hasDemagFFT) {
                        dipolarTemp.x() = dipoleFieldX[site];
                        dipolarTemp.y() = dipoleFieldY[site];
                        dipolarTemp.z() = dipoleFieldZ[site];
                    }

                    effectiveFieldX[site] = exchangeTemp.x() + biasTemp.x() + dipolarTemp.x();
                    effectiveFieldY[site] = exchangeTemp.y() + biasTemp.y() + dipolarTemp.y();
                    effectiveFieldZ[site] = exchangeTemp.z() + biasTemp.z() + dipolarTemp.z();
                }
            }, tbb::auto_partitioner());

    _clearThenReleaseVector(dipoleFieldX);
    _clearThenReleaseVector(dipoleFieldY);
    _clearThenReleaseVector(dipoleFieldZ);
   /*
    if ( simFlags->hasDemagFactors )
        // demagField.DemagnetisationFieldIntense(mxIn, myIn, mzIn, demagXp, demagYp, demagZp);
        demagField.calculateOneDimension(mxIn, myIn, mzIn,
                                         effectiveFieldXAtomic, effectiveFieldYAtomic, effectiveFieldZAtomic,
                                         useParallel);

    _transferDataThenReleaseAtomicVector(effectiveFieldXAtomic, effectiveFieldX);
    _transferDataThenReleaseAtomicVector(effectiveFieldYAtomic, effectiveFieldY);
    _transferDataThenReleaseAtomicVector(effectiveFieldZAtomic, effectiveFieldZ);
    */

    // For data integrity (in case an off-by-one error occurs)
    effectiveFieldX.front() = 0;
    effectiveFieldX.back() = 0;
    effectiveFieldY.front() = 0;
    effectiveFieldY.back() = 0;
    effectiveFieldZ.front() = 0;
    effectiveFieldZ.back() = 0;

    if ( simFlags->hasDipolar )
        dipolarField.DipolarInteractionClassicThreaded(mxIn, myIn, mzIn, dipoleXp, dipoleYp, dipoleZp);

    if ( simFlags->hasSTT )
        // placeholder for example. Find STT for each site at all sites before main loop
        stt.calculateOneDimension(mxIn, myIn, mzIn, sttXp, sttYp, sttZp);

    // Testing. Probably should use single thread, as the overhead here likely won't be worthwhile
    //tbb::global_control c( tbb::global_control::max_allowed_parallelism, 1 );
    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins + 1),
                      [&]( const tbb::blocked_range<int> tbbRange ) {
                          for ( int site = tbbRange.begin(); site < tbbRange.end(); site++ ) {
                              /*
                               * All declarations of overwritten variables are placed here as 'local' to ensure
                               * that each thread has its own definition; else variables are not threadsafe
                              */
                              HkTerms hkLocal;
                              MkTerms mkLocal;

                              // Used for clarity; can be safely refactored away for very minor performance gain
                              hkLocal.x = effectiveFieldX[site];
                              hkLocal.y = effectiveFieldY[site];
                              hkLocal.z = effectiveFieldZ[site];

                              // Calculations for the magnetic moment, coded as symbol 'm', components of the target site
                              auto localGilbertFactor = llg._checkIfDampingMapExists(site);
                              mkLocal.x = llg.MagneticMomentX(site, mxIn[site], myIn[site], mzIn[site],
                                                              hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);
                              mkLocal.y = llg.MagneticMomentY(site, mxIn[site], myIn[site], mzIn[site],
                                                              hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);
                              mkLocal.z = llg.MagneticMomentZ(site, mxIn[site], myIn[site], mzIn[site],
                                                              hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);

                              mxOut[site] = simStates->mx0[site] + stepsize * mkLocal.x;
                              myOut[site] = simStates->my0[site] + stepsize * mkLocal.y;
                              mzOut[site] = simStates->mz0[site] + stepsize * mkLocal.z;

                              _testOutputValues(mxOut[site], myOut[site], mzOut[site], site, iteration, rkStage);
                          }
                      }, tbb::auto_partitioner());

    if (!simFlags->shouldTrackMagneticMomentNorm && rkStage == "2") {
        for (int site = 0; site < simParams->systemTotalSpins; site++) { // Assuming 0-based indexing
            double mIterationNorm = sqrt(pow(mxOut[site], 2) + pow(myOut[site], 2) + pow(mzOut[site], 2));
            double deviation = fabs(1.0 - mIterationNorm); // Absolute deviation from 1.0
            if (deviation > simParams->largestMNorm) {
                simParams->largestMNorm = deviation;
            }
        }
    }

    // if (simParams->largestMNorm > 1.00005) { throw std::runtime_error("mag. moments are no longer below <= 1.00005"); }

    // No need to fill demagX/dipoleX (etc) here they are constantly overwritten; filling to zeroes is just a waste unless debugging

}

void
SolversImplementation::RK4StageMultithreadedCompact( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                                     const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                                     std::vector<double> &myOut, std::vector<double> &mzOut,
                                                     const double &currentTime, const double &stepsize,
                                                     const int &iteration, std::string rkStage ) {
    CommonStructures::Timer dipolarTimer;
    bool useParallel = true;
    int layer = 0;

    // Use these to ensure that there's no issues with sharing a member attribute across
    // methods (this is less memory efficient, but for debugging purposes)
    std::vector<std::atomic<double>> effectiveFieldXAtomic(simParams->systemTotalSpins + 2);
    std::vector<std::atomic<double>> effectiveFieldYAtomic(simParams->systemTotalSpins + 2);
    std::vector<std::atomic<double>> effectiveFieldZAtomic(simParams->systemTotalSpins + 2);

    exchangeField.calculateOneDimension(mxIn, myIn, mzIn, effectiveFieldXAtomic, effectiveFieldYAtomic,
                                        effectiveFieldZAtomic, CommonStructures::Parallelisations::Multithreaded);
    biasField.calculateOneDimension(layer, currentTime, mzIn, effectiveFieldXAtomic, effectiveFieldYAtomic,
                                    effectiveFieldZAtomic, CommonStructures::Parallelisations::Multithreaded);

    if ( simFlags->hasDemagFactors )
        demagField.calculateOneDimension(mxIn, myIn, mzIn, demagXp, demagYp, demagZp,
                                         DemagnetisationFields::DemagMethod::Demag,
                                         CommonStructures::Parallelisations::Multithreaded);

    _transferDataThenReleaseAtomicVector(effectiveFieldXAtomic, effectiveFieldX);
    _transferDataThenReleaseAtomicVector(effectiveFieldYAtomic, effectiveFieldY);
    _transferDataThenReleaseAtomicVector(effectiveFieldZAtomic, effectiveFieldZ);

    if ( simFlags->hasDipolar )
        dipolarField.DipolarInteractionClassicThreaded(mxIn, myIn, mzIn, dipoleXp, dipoleYp, dipoleZp);

    if ( simFlags->hasSTT )
        // placeholder for example. Find STT for each site at all sites before main loop
        stt.calculateOneDimension(mxIn, myIn, mzIn, sttXp, sttYp, sttZp);

    // Testing. Probably should use single thread, as the overhead here likely won't be worthwhile
    //tbb::global_control c( tbb::global_control::max_allowed_parallelism, 1 );
    tbb::parallel_for(tbb::blocked_range<int>(1, simParams->systemTotalSpins + 1),
                      [&]( const tbb::blocked_range<int> tbbRange ) {
                          for ( int site = tbbRange.begin(); site < tbbRange.end(); site++ ) {
                              /*
                               * All declarations of overwritten variables are placed here as 'local' to ensure
                               * that each thread has its own definition; else variables are not threadsafe
                              */
                              HkTerms hkLocal;
                              MkTerms mkLocal;

                              // Used for clarity; can be safely refactored away for very minor performance gain
                              hkLocal.x = effectiveFieldX[site];
                              hkLocal.y = effectiveFieldY[site];
                              hkLocal.z = effectiveFieldZ[site];

                              // Calculations for the magnetic moment, coded as symbol 'm', components of the target site
                              double localGilbertFactor = llg._checkIfDampingMapExists(site);
                              mkLocal.x = llg.MagneticMomentX(site, mxIn[site], myIn[site], mzIn[site],
                                                              hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);
                              mkLocal.y = llg.MagneticMomentY(site, mxIn[site], myIn[site], mzIn[site],
                                                              hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);
                              mkLocal.z = llg.MagneticMomentZ(site, mxIn[site], myIn[site], mzIn[site],
                                                              hkLocal.x, hkLocal.y, hkLocal.z, localGilbertFactor);

                              if ( rkStage == "2" or rkStage == "3" ) {
                                  mxk[site] += 2 * mkLocal.x;
                                  myk[site] += 2 * mkLocal.y;
                                  mzk[site] += 2 * mkLocal.z;
                              } else {
                                  mxk[site] += mkLocal.x;
                                  myk[site] += mkLocal.y;
                                  mzk[site] += mkLocal.z;
                              }

                              mxOut[site] = simStates->mx0[site] + stepsize * mkLocal.x;
                              myOut[site] = simStates->my0[site] + stepsize * mkLocal.y;
                              mzOut[site] = simStates->mz0[site] + stepsize * mkLocal.z;

                              _testOutputValues(mxOut[site], myOut[site], mzOut[site], site, iteration, rkStage);
                          }
                      }, tbb::auto_partitioner());

    if (!simFlags->shouldTrackMagneticMomentNorm && rkStage == "2") {
        for (int site = 0; site < simParams->systemTotalSpins; site++) { // Assuming 0-based indexing
            double mIterationNorm = sqrt(pow(mxOut[site], 2) + pow(myOut[site], 2) + pow(mzOut[site], 2));
            double deviation = fabs(1.0 - mIterationNorm); // Absolute deviation from 1.0
            if (deviation > simParams->largestMNorm) {
                simParams->largestMNorm = deviation;
            }
        }
    }

    // if (simParams->largestMNorm > 1.00005) { throw std::runtime_error("mag. moments are no longer below <= 1.00005"); }

    // No need to fill demagX/dipoleX (etc) here they are constantly overwritten; filling to zeroes is just a waste unless debugging

}

void SolversImplementation::runMethod() {

    std::string methodToUse = GV.GetNumericalMethod();

    if ( methodToUse == "RK2" )
        SolveRK2();
    else if ( methodToUse == "RK2c" )
        SolveRK2Classic();
    else if ( methodToUse == "RK2p" )
        RK2Parallel();
    else if ( methodToUse == "RK4p" )
        RK4Parallel();
    else
        throw std::runtime_error("Method not recognised");
}

void SolversImplementation::_resizeClassContainersTest() {
    // Should contain all interactions/fields that are calculated
    effectiveFieldX.resize(simParams->systemTotalSpins + 2);
    effectiveFieldY.resize(simParams->systemTotalSpins + 2);
    effectiveFieldZ.resize(simParams->systemTotalSpins + 2);

    mx1p.resize(simParams->systemTotalSpins + 2);
    my1p.resize(simParams->systemTotalSpins + 2);
    mz1p.resize(simParams->systemTotalSpins + 2);
    mx2p.resize(simParams->systemTotalSpins + 2);
    my2p.resize(simParams->systemTotalSpins + 2);
    mz2p.resize(simParams->systemTotalSpins + 2);
}

void SolversImplementation::_resizeClassContainers() {
    // Should contain all interactions/fields that are calculated

    if ((simFlags->hasDemagFactors) || (!simFlags->hasDemagFactors)) {
        demagXp.resize(simParams->systemTotalSpins + 2);
        std::fill(demagXp.begin(), demagXp.end(), 0.0);
        demagYp.resize(simParams->systemTotalSpins + 2);
        std::fill(demagYp.begin(), demagYp.end(), 0.0);
        demagZp.resize(simParams->systemTotalSpins + 2);
        std::fill(demagZp.begin(), demagZp.end(), 0.0);
    }
    // TODO. Fix ugly, nasty code below
    // Horrible IF statement is because not all methods currently access single elements and instead use the whole vector
    // There would then be a SIGSEGV error if the vector was not resized even if the vector was never used during a calculation
    if ( simFlags->hasDipolar ) {
        dipoleXp.resize(simParams->systemTotalSpins + 2);
        std::fill(dipoleXp.begin(), dipoleXp.end(), 0.0);
        dipoleYp.resize(simParams->systemTotalSpins + 2);
        std::fill(dipoleYp.begin(), dipoleYp.end(), 0.0);
        dipoleZp.resize(simParams->systemTotalSpins + 2);
        std::fill(dipoleZp.begin(), dipoleZp.end(), 0.0);
    }

    if ( simFlags->hasDMI ) {
        dmiXp.resize(simParams->systemTotalSpins + 2);
        std::fill(dmiXp.begin(), dmiXp.end(), 0.0);
        dmiYp.resize(simParams->systemTotalSpins + 2);
        std::fill(dmiYp.begin(), dmiYp.end(), 0.0);
        dmiZp.resize(simParams->systemTotalSpins + 2);
        std::fill(dmiZp.begin(), dmiZp.end(), 0.0);
    }

    if ( simFlags->hasSTT ) {
        sttXp.resize(simParams->systemTotalSpins + 2);
        std::fill(sttXp.begin(), sttXp.end(), 0.0);
        sttYp.resize(simParams->systemTotalSpins + 2);
        std::fill(sttYp.begin(), sttYp.end(), 0.0);
        sttZp.resize(simParams->systemTotalSpins + 2);
        std::fill(sttZp.begin(), sttZp.end(), 0.0);
    }

    // Fill RK2 Stage magnetic moment containers. Ensure only to include methods that use class-wide containers
    if ( GV.GetNumericalMethod() == "RK2p" ) {
        mx1p.resize(simParams->systemTotalSpins + 2);
        my1p.resize(simParams->systemTotalSpins + 2);
        mz1p.resize(simParams->systemTotalSpins + 2);
        mx2p.resize(simParams->systemTotalSpins + 2);
        my2p.resize(simParams->systemTotalSpins + 2);
        mz2p.resize(simParams->systemTotalSpins + 2);
    }
}

void SolversImplementation::_resizeClassContainersRK4() {
    // Should contain all interactions/fields that are calculated
    effectiveFieldX.resize(simParams->systemTotalSpins + 2);
    effectiveFieldY.resize(simParams->systemTotalSpins + 2);
    effectiveFieldZ.resize(simParams->systemTotalSpins + 2);

    mx1p.resize(simParams->systemTotalSpins + 2);
    my1p.resize(simParams->systemTotalSpins + 2);
    mz1p.resize(simParams->systemTotalSpins + 2);
    mx2p.resize(simParams->systemTotalSpins + 2);
    my2p.resize(simParams->systemTotalSpins + 2);
    mz2p.resize(simParams->systemTotalSpins + 2);
    mx3p.resize(simParams->systemTotalSpins + 2);
    my3p.resize(simParams->systemTotalSpins + 2);
    mz3p.resize(simParams->systemTotalSpins + 2);
    mx4p.resize(simParams->systemTotalSpins + 2);
    my4p.resize(simParams->systemTotalSpins + 2);
    mz4p.resize(simParams->systemTotalSpins + 2);

    mxk.resize(simParams->systemTotalSpins + 2);
    myk.resize(simParams->systemTotalSpins + 2);
    mzk.resize(simParams->systemTotalSpins + 2);
}

void SolversImplementation::_clearDataStorageVariables() {

    simParams->totalTime = 0.0;
    std::fill(mx1p.begin(), mx1p.end(), 0.0);
    std::fill(my1p.begin(), my1p.end(), 0.0);
    std::fill(mz1p.begin(), mz1p.end(), 0.0);
    std::fill(mx2p.begin(), mx2p.end(), 0.0);
    std::fill(my2p.begin(), my2p.end(), 0.0);
    std::fill(mz2p.begin(), mz2p.end(), 0.0);
}

void SolversImplementation::_testOutputValues( double &mxTerm, double &myTerm, double &mzTerm, int site, int iteration,
                                               const std::string &rkStage ) {
    if ( std::isinf(mxTerm))
        throw std::runtime_error(
                "mxOut is INF at site " + std::to_string(site) + " at iteration " +
                std::to_string(iteration) + " in RK2 stage " + rkStage);
    if ( std::isnan(mxTerm))
        throw std::runtime_error(
                "mxOut is NaN at site " + std::to_string(site) + " at iteration " +
                std::to_string(iteration) + " in RK2 stage " + rkStage);

    if ( std::isinf(myTerm))
        throw std::runtime_error(
                "myOut is INF at site " + std::to_string(site) + " at iteration " +
                std::to_string(iteration) + " in RK2 stage " + rkStage);
    if ( std::isnan(myTerm))
        throw std::runtime_error(
                "myOut is NaN at site " + std::to_string(site) + " at iteration " +
                std::to_string(iteration) + " in RK2 stage " + rkStage);

    if ( std::isinf(mzTerm))
        throw std::runtime_error(
                "mzOut is INF at site " + std::to_string(site) + " at iteration " +
                std::to_string(iteration) + " in RK2 stage " + rkStage);
    if ( std::isnan(mzTerm))
        throw std::runtime_error(
                "mzOut is NaN at site " + std::to_string(site) + " at iteration " +
                std::to_string(iteration) + " in RK2 stage " + rkStage);
}

void SolversImplementation::_transferDataThenReleaseAtomicVector( std::vector<std::atomic<double>> &atomicVector,
                                                                  std::vector<double> &regularVector,
                                                                  bool shouldRelease ) {
    for ( size_t i = 0; i < atomicVector.size(); i++ )
        regularVector[i] = atomicVector[i].load();

    if ( shouldRelease ) {
        atomicVector.clear();
        atomicVector.shrink_to_fit();
    }

    /*
    for (size_t i = 1; i <= simParams->systemTotalSpins; i++) {
        effectiveFieldX[i] = effectiveFieldXLocal[i].load();
        effectiveFieldY[i] = effectiveFieldYLocal[i].load();
        effectiveFieldZ[i] = effectiveFieldZLocal[i].load();
    }
     */
}

void SolversImplementation::_clearThenReleaseVector( auto &vec ) {
    vec.clear();
    vec.shrink_to_fit();
}