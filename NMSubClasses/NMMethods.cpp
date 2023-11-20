//
// Created by Cameron McEleney on 31/10/2023.
//

// C++ User Libraries (Current)
#include "NMMethods.h"


// C++ User Libraries (Children)

NMMethods::NMMethods(std::shared_ptr<SystemDataContainer> data) :
        NMSuperClassTest(std::move(data)), 
        demagField(simState.get()),
        effectiveField(simState.get()),
        dipolarField(simState.get()),
        llg(simState.get()){}

void NMMethods::_testShockwaveConditions(double iteration) {

    if (simState->shouldDriveCease) {
        // and (simState->isShockwaveOn and simState->isShockwaveAtMax)) {
        if (simState->isShockwaveOn and not simState->isShockwaveAtMax)
            std::cout << "Shock not at maximum when cut-off" << std::endl;

        if (iteration >= simState->iterationEnd * simState->iterEndShock) {
            // Shockwave begins once simulation is a certain % complete
            simState->hasShockwave = false;
            simState->isShockwaveOn = false;
            simState->dynamicBiasField = 0;
        }

        return;
    }

    // If method is triggered, then the applied biasFieldDriving is increased by the scale factor shockwaveScaling
    if (simState->hasShockwave and not simState->isShockwaveOn)
    {
        if (iteration >= simState->iterationEnd * simState->iterStartShock)
        {
            // Shockwave begins once simulation is a certain % complete
            simState->isShockwaveOn = true;
            simState->dynamicBiasField = simState->shockwaveInitialStrength;
        }

        return;
    }

    if (simState->isShockwaveOn and not simState->isShockwaveAtMax)
    {
        simState->dynamicBiasField += simState->shockwaveStepsize;

        if (simState->dynamicBiasField >= simState->shockwaveMax)
        {
            simState->dynamicBiasField = simState->shockwaveMax;
            simState->isShockwaveAtMax = true;

        }
        return;

    }

}

void NMMethods::SolveRK2Classic() {
    std::shared_ptr<NMDataHandling> childNMData = std::make_shared<NMDataHandling>(simState);
    // Only uses a single spin chain to solve the RK2 midpoint method.

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    //std::ofstream myRK2File(GV.GetFilePath() + "rk2_my_" + GV.GetFileNameBase() + ".csv");
    //std::ofstream mzRK2File(GV.GetFilePath() + "rk2_mz_" + GV.GetFileNameBase() + ".csv");

    std::cout << "test val of drivingFreq: " << simState->drivingFreq << std::endl;
    std::cout << "Trying to use childNMData" << std::endl;
    if (simState->isFm) {
        childNMData->InformUserOfCodeType("RK2 Midpoint (FM)");
        childNMData->CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)");
        //childNMData -> CreateFileHeader(myRK2File, "RK2 Midpoint (FM)");
        //childNMData -> CreateFileHeader(mzRK2File, "RK2 Midpoint (FM)");
    std::cout << "Managed to use childNMData" << std::endl;
    } else if (!simState->isFm) {
        childNMData -> InformUserOfCodeType("RK2 Midpoint (AFM)");
        childNMData -> CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
    }

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata();
    }

    progressbar bar(100);

    std::vector<double> demagX(GV.GetNumSpins() + 2, 0.0), demagY(GV.GetNumSpins() + 2, 0.0), demagZ(GV.GetNumSpins() + 2, 0.0);
    std::vector<double> dipoleX(GV.GetNumSpins() + 2, 0.0), dipoleY(GV.GetNumSpins() + 2, 0.0), dipoleZ(GV.GetNumSpins() + 2, 0.0);

    for (int iteration = simState->iterationStart; iteration <= simState->iterationEnd; iteration++) {

        if (simState->iterationEnd >= 100 && iteration % (simState->iterationEnd / 100) == 0)
            bar.update(); // Doesn't work for fewer than 100 iterations

        _testShockwaveConditions(iteration);

        double t0 = simState->totalTime;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = simState->mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(GV.GetNumSpins() + 2, 0), my1(GV.GetNumSpins() + 2, 0), mz1(GV.GetNumSpins() + 2, 0);
        // EASY FIND
        if (simState->useDipolar)
            dipolarField.DipolarInteraction1D(simState->mx0, dipoleX);
        if (simState->useDemagIntense) {
            demagField.DemagnetisationFieldIntense(demagX, demagY, demagZ, simState->mx0, simState->my0, simState->mz0);
        } else if (simState->useDemagFft) {
            std::string rkStageName = "2-1";
            demagField.DemagField1DReal(demagX, demagY, demagZ, simState->mx0, simState->my0, simState->mz0, iteration, rkStageName);

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

        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;
            /*
            if ((simState->iterationEnd >= 100 && iteration % (simState->iterationEnd / 1000) == 0) && site == 500) {
                std::cout << "Iter. #" << iteration << " ";
                std::cout << "| mx: " << simState->mx0[200] << " - H_dx: " << demagX[200] << " ";
                std::cout << "| my: " << simState->my0[200] << " - H_dy: " << demagY[200] << " ";
                std::cout << "| mz: " << simState->mz0[200] << " - H_dz: " << demagZ[200] << " ";
                std::cout << "| mTot: " << sqrt(pow(simState->mx0[site], 2) + pow(simState->my0[site], 2) + pow(simState->mz0[site], 2)) << " ";
                std::cout << std::endl;
                //std::cout << "| mz: " << simState->mz0[200] << " ";
                //std::cout << "| H_dz: " << demagZ[200] << " ";
                //std::cout << "| mTot: " << sqrt(pow(simState->mx0[site], 2) + pow(simState->my0[site], 2) + pow(simState->mz0[site], 2)) << " ";
                //std::cout << std::endl;
            }
            */
            /*double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if (simState->useDipolar) {
                std::vector<double> mxTermsForDipole = {simState->mx0[spinLHS], simState->mx0[site], simState->mx0[spinRHS]};
                std::vector<double> myTermsForDipole = {simState->my0[spinLHS], simState->my0[site], simState->my0[spinRHS]};
                std::vector<double> mzTermsForDipole = {simState->mz0[spinLHS], simState->mz0[site], simState->mz0[spinRHS]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipolarInteractionClassic(mxTermsForDipole, myTermsForDipole,
                                                                            mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            }
            */
            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK0 = effectiveField.EffectiveFieldX(site, 0, simState->mx0[spinLHS], simState->mx0[site], simState->mx0[spinRHS], dipoleX[site], demagX[site], t0);
            double hyK0 = effectiveField.EffectiveFieldY(site, 0, simState->my0[spinLHS], simState->my0[site], simState->my0[spinRHS], dipoleY[site], demagY[site]);
            double hzK0 = effectiveField.EffectiveFieldZ(site, 0, simState->mz0[spinLHS], simState->mz0[site], simState->mz0[spinRHS], dipoleZ[site], demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK1 = llg.MagneticMomentX(site, simState->mx0[site], simState->my0[site], simState->mz0[site], hxK0, hyK0, hzK0);
            double myK1 = llg.MagneticMomentY(site, simState->mx0[site], simState->my0[site], simState->mz0[site], hxK0, hyK0, hzK0);
            double mzK1 = llg.MagneticMomentZ(site, simState->mx0[site], simState->my0[site], simState->mz0[site], hxK0, hyK0, hzK0);

            // Find (m0 + k1/2) for each site, which is used in the next stage.
            mx1[site] = simState->mx0[site] + simState->stepsizeHalf * mxK1;
            my1[site] = simState->my0[site] + simState->stepsizeHalf * myK1;
            mz1[site] = simState->mz0[site] + simState->stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        // EASY FIND
        std::vector<double> mx2(GV.GetNumSpins() + 2, 0), my2(GV.GetNumSpins() + 2, 0), mz2(GV.GetNumSpins() + 2, 0);
        std::fill(demagX.begin(), demagX.end(), 0.0); std::fill(demagY.begin(), demagY.end(), 0.0); std::fill(demagZ.begin(), demagZ.end(), 0.0);
        std::fill(dipoleX.begin(), dipoleX.end(), 0.0); std::fill(dipoleY.begin(), dipoleY.end(), 0.0); std::fill(dipoleZ.begin(), dipoleZ.end(), 0.0);

        if (simState->useDipolar)
            dipolarField.DipolarInteraction1D(mx1, dipoleX);
        if (simState->useDemagIntense) {
            demagField.DemagnetisationFieldIntense(demagX, demagY, demagZ, mx1, my1, mz1);
        } else if (simState->useDemagFft) {
            std::string rkStageName = "2-2";
            demagField.DemagField1DReal(demagX, demagY, demagZ, mx1, my1, mz1, iteration, rkStageName);
            // if (iteration > 0) {std::cout << "Stage 2" << std::endl; PrintVector(demagZ, false);}
        }

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            /*double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if (simState->useDipolar) {
                std::vector<double> mxTermsForDipole = {mx1[spinLHS], mx1[site], mx1[spinRHS]};
                std::vector<double> myTermsForDipole = {my1[spinLHS], my1[site], my1[spinRHS]};
                std::vector<double> mzTermsForDipole = {mz1[spinLHS], mz1[site], mz1[spinRHS]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipolarInteractionClassic(mxTermsForDipole, myTermsForDipole,
                                                                       mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            }
            */
            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK1 = effectiveField.EffectiveFieldX(site, 0, mx1[spinLHS], mx1[site], mx1[spinRHS], dipoleX[site], demagX[site], t0);
            double hyK1 = effectiveField.EffectiveFieldY(site, 0, my1[spinLHS], my1[site], my1[spinRHS], dipoleY[site], demagY[site]);
            double hzK1 = effectiveField.EffectiveFieldZ(site, 0, mz1[spinLHS], mz1[site], mz1[spinRHS], dipoleZ[site], demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK2 = llg.MagneticMomentX(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double myK2 = llg.MagneticMomentY(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double mzK2 = llg.MagneticMomentZ(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);

            mx2[site] = simState->mx0[site] + simState->stepsize * mxK2;
            my2[site] = simState->my0[site] + simState->stepsize * myK2;
            mz2[site] = simState->mz0[site] + simState->stepsize * mzK2;

            if (simState->shouldTrackMValues) {
                double mIterationNorm = sqrt(pow(mx2[site], 2) + pow(my2[site], 2) + pow(mz2[site], 2));
                if ((simState->largestMNorm) > (1.0 - mIterationNorm)) { simState->largestMNorm = (1.0 - mIterationNorm); }
                if (mIterationNorm > 1.00005) {throw std::runtime_error("mag. moments are no longer below <= 1.00005");}
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        simState->mx0.clear(); simState->my0.clear(); simState->mz0.clear();
        mx1.clear(); my1.clear(); mz1.clear();

        childNMData -> SaveDataToFile(mxRK2File, mx2, iteration);
        //childNMData -> SaveDataToFile(myRK2File, my2, iteration);
        //childNMData -> SaveDataToFile(mzRK2File, mz2, iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        simState->mx0 = mx2; simState->my0 = my2; simState->mz0 = mz2;

        if (iteration == simState->forceStopAtIteration)
            exit(0);

        simState->totalTime += simState->stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();
    //myRK2File.close();
    //mzRK2File.close();

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata(true);
    }

    if (simState->shouldTrackMValues)
        std::cout << "\nMax norm. value of M is: " << simState->largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}
void NMMethods::SolveRK2() {
    std::shared_ptr<NMDataHandling> childNMData = std::make_shared<NMDataHandling>(simState);
    // Uses multiple layers to solve the RK2 midpoint method. See the documentation for more details.

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    std::ofstream mxRK2File1(GV.GetFilePath() + "rk2_mx1_" + GV.GetFileNameBase() + ".csv");

    // User information and file header is magnetic-material specific.
    if (simState->isFm) {
        childNMData -> InformUserOfCodeType("RK2 Midpoint (FM)");
        childNMData -> CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)", false, 0);
        childNMData -> CreateFileHeader(mxRK2File1, "RK2 Midpoint (FM)", false, 1);
    } else if (!simState->isFm) {
        childNMData -> InformUserOfCodeType("RK2 Midpoint (AFM)");
        childNMData -> CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
        childNMData -> CreateFileHeader(mxRK2File1, "RK2 Midpoint (AFM)");
    }

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata();
    }

    progressbar bar(100);

    std::vector<double> demagX(GV.GetNumSpins() + 2, 0.0), demagY(GV.GetNumSpins() + 2, 0.0), demagZ(GV.GetNumSpins() + 2, 0.0);

    for (int iteration = simState->iterationStart; iteration <= simState->iterationEnd; iteration++) {

        if (simState->iterationEnd >= 100 && iteration % (simState->iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        _testShockwaveConditions(iteration);

        double t0 = simState->totalTime;

        for (int layer = 0; layer < simState->totalLayers; layer++) {
            // RK2 Stage 1. Takes initial conditions as inputs.

            for (int site = 1; site <= simState->layerTotalSpins[layer]; site++) {
                // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)

                // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
                int spinLHS = site - 1, spinRHS = site + 1;

                double mxLHS = simState->m0Nest[layer][spinLHS][0], mxMID = simState->m0Nest[layer][site][0], mxRHS = simState->m0Nest[layer][spinRHS][0];
                double myLHS = simState->m0Nest[layer][spinLHS][1], myMID = simState->m0Nest[layer][site][1], myRHS = simState->m0Nest[layer][spinRHS][1];
                double mzLHS = simState->m0Nest[layer][spinLHS][2], mzMID = simState->m0Nest[layer][site][2], mzRHS = simState->m0Nest[layer][spinRHS][2];

                double dipoleX, dipoleY, dipoleZ;
                if (simState->useDipolar) {

                    int layer1, layer2;
                    if (layer == 0) {layer1 = 0; layer2 = 1;}
                    else if (layer == 1) {layer1 = 1; layer2 = 0;}

                    if (simState->debugFunc) {std::cout << "\n\niteration: " << iteration << " | layer: " << layer << " | site: " << site << std::endl;}
                    std::vector<double> dipoleTerms = dipolarField.DipolarInteractionInterlayer(simState->m0Nest[layer1], simState->m0Nest[layer2], site,
                                                                                                layer1, layer2);

                    dipoleX = dipoleTerms[0];
                    dipoleY = dipoleTerms[1];
                    dipoleZ = dipoleTerms[2];
                } else {
                    dipoleX = 0;
                    dipoleY = 0;
                    dipoleZ = 0;
                }

                // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                double hxK0 = effectiveField.EffectiveFieldX(site, layer, mxLHS, mxMID, mxRHS, dipoleX, demagX[site], t0);
                double hyK0 = effectiveField.EffectiveFieldY(site, layer, myLHS, myMID, myRHS, dipoleY, demagY[site]);
                double hzK0 = effectiveField.EffectiveFieldZ(site, layer, mzLHS, mzMID, mzRHS, dipoleZ, demagZ[site]);

                // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
                double mxK1 = llg.MagneticMomentX(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
                double myK1 = llg.MagneticMomentY(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
                double mzK1 = llg.MagneticMomentZ(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);

                // Find (m0 + k1/2) for each site, which is used in the next stage.
                simState->m1Nest[layer][site][0] = mxMID + simState->stepsizeHalf * mxK1;
                simState->m1Nest[layer][site][1] = myMID + simState->stepsizeHalf * myK1;
                simState->m1Nest[layer][site][2] = mzMID + simState->stepsizeHalf * mzK1;
            }
        }

        for (int layer = 0; layer < simState->totalLayers; layer++) {
            // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
            for (int site = 1; site <= simState->layerTotalSpins[layer]; site++) {

                // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
                int spinLHS = site - 1, spinRHS = site + 1;

                double mxLHS = simState->m1Nest[layer][spinLHS][0], mxMID = simState->m1Nest[layer][site][0], mxRHS = simState->m1Nest[layer][spinRHS][0];
                double myLHS = simState->m1Nest[layer][spinLHS][1], myMID = simState->m1Nest[layer][site][1], myRHS = simState->m1Nest[layer][spinRHS][1];
                double mzLHS = simState->m1Nest[layer][spinLHS][2], mzMID = simState->m1Nest[layer][site][2], mzRHS = simState->m1Nest[layer][spinRHS][2];

                double dipoleX, dipoleY, dipoleZ;
                if (simState->useDipolar) {

                    int layer1, layer2;
                    if (layer == 0) {layer1 = 0; layer2 = 1;}
                    else if (layer == 1) {layer1 = 1; layer2 = 0;}

                    int debugCounter = 0;  // To make sure debug outputs only occur during the first RK2 stage, not this second stage
                    if (simState->debugFunc) { simState->debugFunc = false; debugCounter++; }
                    std::vector<double> dipoleTerms = dipolarField.DipolarInteractionInterlayer(simState->m1Nest[layer1], simState->m1Nest[layer2], site,
                                                                                                layer1, layer2);
                    if (debugCounter > 0) { simState->debugFunc = true; }

                    dipoleX = dipoleTerms[0];
                    dipoleY = dipoleTerms[1];
                    dipoleZ = dipoleTerms[2];
                } else {
                    dipoleX = 0;
                    dipoleY = 0;
                    dipoleZ = 0;
                }
                // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                double hxK1 = effectiveField.EffectiveFieldX(site, layer, mxLHS, mxMID, mxRHS, dipoleX, demagX[site], t0);
                double hyK1 = effectiveField.EffectiveFieldY(site, layer, myLHS, myMID, myRHS, dipoleY, demagY[site]);
                double hzK1 = effectiveField.EffectiveFieldZ(site, layer, mzLHS, mzMID, mzRHS, dipoleZ, demagZ[site]);

                // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
                double mxK2 = llg.MagneticMomentX(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
                double myK2 = llg.MagneticMomentY(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
                double mzK2 = llg.MagneticMomentZ(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);

                simState->m2Nest[layer][site][0] = simState->m0Nest[layer][site][0] + simState->stepsize * mxK2;
                simState->m2Nest[layer][site][1] = simState->m0Nest[layer][site][1] + simState->stepsize * myK2;
                simState->m2Nest[layer][site][2] = simState->m0Nest[layer][site][2] + simState->stepsize * mzK2;

                if (simState->shouldTrackMValues) {
                    double mIterationNorm = sqrt(
                            pow(simState->m2Nest[layer][site][0], 2) + pow(simState->m2Nest[layer][site][1], 2) + pow(simState->m2Nest[layer][site][2], 2));
                    if ((simState->largestMNormMulti[layer]) > (1.0 - mIterationNorm)) { simState->largestMNormMulti[layer] = (1.0 - mIterationNorm); }
                }
            }
        }

        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */

        childNMData -> SaveDataToFileMultilayer(mxRK2File, simState->m2Nest[0], iteration, 0);
        childNMData -> SaveDataToFileMultilayer(mxRK2File1, simState->m2Nest[1], iteration, 1);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        simState->m0Nest = simState->m2Nest;

        if (iteration == simState->forceStopAtIteration)
            exit(0);

        simState->totalTime += simState->stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata(true);
    }

    if (simState->shouldTrackMValues) {
        std::cout << "\nMax norm. values of M are: ";
        for (int i = 0; i < simState->largestMNormMulti.size(); i++) {
            std::cout << "Layer " << i << ": " << simState->largestMNormMulti[i] << " | ";
        }
    }

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}

void NMMethods::runMethod() {

    std::cout<<"Running method" << std::endl;

    std::string methodToUse = "RK2c";
    if (methodToUse == "RK2")
        SolveRK2();
    else if (methodToUse == "RK2c")
        SolveRK2Classic();
    else
        throw std::runtime_error("Method not recognised");
}