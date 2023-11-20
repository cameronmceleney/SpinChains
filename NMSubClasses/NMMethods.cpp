//
// Created by Cameron McEleney on 31/10/2023.
//

// C++ User Libraries (Current)
#include "NMMethods.h"

#include <utility>

#include "SystemDataContainer.h"

// C++ User Libraries (Sibling Classes)
#include "NMDataHandling.h"

// C++ User Libraries (Children)

NMMethods::NMMethods(std::shared_ptr<SystemDataContainer> data) :
        NMSuperClassTest(std::move(data))
        {}

void NMMethods::_testShockwaveConditions(double iteration) {

    if (systemData->shouldDriveCease) {
        // and (systemData->isShockwaveOn and systemData->isShockwaveAtMax)) {
        if (systemData->isShockwaveOn and not systemData->isShockwaveAtMax)
            std::cout << "Shock not at maximum when cut-off" << std::endl;

        if (iteration >= systemData->iterationEnd * systemData->iterEndShock) {
            // Shockwave begins once simulation is a certain % complete
            systemData->hasShockwave = false;
            systemData->isShockwaveOn = false;
            systemData->dynamicBiasField = 0;
        }

        return;
    }

    // If method is triggered, then the applied biasFieldDriving is increased by the scale factor shockwaveScaling
    if (systemData->hasShockwave and not systemData->isShockwaveOn)
    {
        if (iteration >= systemData->iterationEnd * systemData->iterStartShock)
        {
            // Shockwave begins once simulation is a certain % complete
            systemData->isShockwaveOn = true;
            systemData->dynamicBiasField = systemData->shockwaveInitialStrength;
        }

        return;
    }

    if (systemData->isShockwaveOn and not systemData->isShockwaveAtMax)
    {
        systemData->dynamicBiasField += systemData->shockwaveStepsize;

        if (systemData->dynamicBiasField >= systemData->shockwaveMax)
        {
            systemData->dynamicBiasField = systemData->shockwaveMax;
            systemData->isShockwaveAtMax = true;

        }
        return;

    }

}

void NMMethods::SolveRK2Classic() {
    std::shared_ptr<NMDataHandling> childNMData = std::make_shared<NMDataHandling>(systemData);
    // Only uses a single spin chain to solve the RK2 midpoint method.

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    //std::ofstream myRK2File(GV.GetFilePath() + "rk2_my_" + GV.GetFileNameBase() + ".csv");
    //std::ofstream mzRK2File(GV.GetFilePath() + "rk2_mz_" + GV.GetFileNameBase() + ".csv");

    std::cout << "test val of drivingFreq: " << systemData->drivingFreq << std::endl;
    std::cout << "Trying to use childNMData" << std::endl;
    if (systemData->isFm) {
        childNMData->InformUserOfCodeType("RK2 Midpoint (FM)");
        childNMData->CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)");
        //childNMData -> CreateFileHeader(myRK2File, "RK2 Midpoint (FM)");
        //childNMData -> CreateFileHeader(mzRK2File, "RK2 Midpoint (FM)");
    std::cout << "Managed to use childNMData" << std::endl;

    } else if (!systemData->isFm) {
        childNMData -> InformUserOfCodeType("RK2 Midpoint (AFM)");
        childNMData -> CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
    }

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata();
    }

    progressbar bar(100);

    std::vector<double> demagX(GV.GetNumSpins() + 2, 0.0), demagY(GV.GetNumSpins() + 2, 0.0), demagZ(GV.GetNumSpins() + 2, 0.0);
    std::vector<double> dipoleX(GV.GetNumSpins() + 2, 0.0), dipoleY(GV.GetNumSpins() + 2, 0.0), dipoleZ(GV.GetNumSpins() + 2, 0.0);

    for (int iteration = systemData->iterationStart; iteration <= systemData->iterationEnd; iteration++) {

        if (systemData->iterationEnd >= 100 && iteration % (systemData->iterationEnd / 100) == 0)
            bar.update(); // Doesn't work for fewer than 100 iterations

        _testShockwaveConditions(iteration);

        double t0 = systemData->totalTime;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = systemData->mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(GV.GetNumSpins() + 2, 0), my1(GV.GetNumSpins() + 2, 0), mz1(GV.GetNumSpins() + 2, 0);
        // EASY FIND
        if (systemData->useDipolar)
            Dipolar -> DipolarInteraction1D(systemData->mx0, dipoleX);
        if (systemData->useDemagIntense) {
            demagField -> DemagnetisationFieldIntense(demagX, demagY, demagZ, systemData->mx0, systemData->my0, systemData->mz0);
        } else if (systemData->useDemagFft) {
            std::string rkStageName = "2-1";
            demagField -> DemagField1DReal(demagX, demagY, demagZ, systemData->mx0, systemData->my0, systemData->mz0, iteration, rkStageName);

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
            if ((systemData->iterationEnd >= 100 && iteration % (systemData->iterationEnd / 1000) == 0) && site == 500) {
                std::cout << "Iter. #" << iteration << " ";
                std::cout << "| mx: " << systemData->mx0[200] << " - H_dx: " << demagX[200] << " ";
                std::cout << "| my: " << systemData->my0[200] << " - H_dy: " << demagY[200] << " ";
                std::cout << "| mz: " << systemData->mz0[200] << " - H_dz: " << demagZ[200] << " ";
                std::cout << "| mTot: " << sqrt(pow(systemData->mx0[site], 2) + pow(systemData->my0[site], 2) + pow(systemData->mz0[site], 2)) << " ";
                std::cout << std::endl;
                //std::cout << "| mz: " << systemData->mz0[200] << " ";
                //std::cout << "| H_dz: " << demagZ[200] << " ";
                //std::cout << "| mTot: " << sqrt(pow(systemData->mx0[site], 2) + pow(systemData->my0[site], 2) + pow(systemData->mz0[site], 2)) << " ";
                //std::cout << std::endl;
            }
            */
            /*double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if (systemData->useDipolar) {
                std::vector<double> mxTermsForDipole = {systemData->mx0[spinLHS], systemData->mx0[site], systemData->mx0[spinRHS]};
                std::vector<double> myTermsForDipole = {systemData->my0[spinLHS], systemData->my0[site], systemData->my0[spinRHS]};
                std::vector<double> mzTermsForDipole = {systemData->mz0[spinLHS], systemData->mz0[site], systemData->mz0[spinRHS]};
                std::vector<int> siteTermsForDipole = {spinLHS, site, spinRHS};

                std::vector<double> dipoleTerms = DipolarInteractionClassic(mxTermsForDipole, myTermsForDipole,
                                                                            mzTermsForDipole, siteTermsForDipole);

                dipoleX = dipoleTerms[0];
                dipoleY = dipoleTerms[1];
                dipoleZ = dipoleTerms[2];
            }
            */
            // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
            double hxK0 = EffField->EffectiveFieldX(site, 0, systemData->mx0[spinLHS], systemData->mx0[site], systemData->mx0[spinRHS], dipoleX[site], demagX[site], t0);
            double hyK0 = EffField -> EffectiveFieldY(site, 0, systemData->my0[spinLHS], systemData->my0[site], systemData->my0[spinRHS], dipoleY[site], demagY[site]);
            double hzK0 = EffField -> EffectiveFieldZ(site, 0, systemData->mz0[spinLHS], systemData->mz0[site], systemData->mz0[spinRHS], dipoleZ[site], demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK1 = LLG -> MagneticMomentX(site, systemData->mx0[site], systemData->my0[site], systemData->mz0[site], hxK0, hyK0, hzK0);
            double myK1 = LLG -> MagneticMomentY(site, systemData->mx0[site], systemData->my0[site], systemData->mz0[site], hxK0, hyK0, hzK0);
            double mzK1 = LLG -> MagneticMomentZ(site, systemData->mx0[site], systemData->my0[site], systemData->mz0[site], hxK0, hyK0, hzK0);

            // Find (m0 + k1/2) for each site, which is used in the next stage.
            mx1[site] = systemData->mx0[site] + systemData->stepsizeHalf * mxK1;
            my1[site] = systemData->my0[site] + systemData->stepsizeHalf * myK1;
            mz1[site] = systemData->mz0[site] + systemData->stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        // EASY FIND
        std::vector<double> mx2(GV.GetNumSpins() + 2, 0), my2(GV.GetNumSpins() + 2, 0), mz2(GV.GetNumSpins() + 2, 0);
        std::fill(demagX.begin(), demagX.end(), 0.0); std::fill(demagY.begin(), demagY.end(), 0.0); std::fill(demagZ.begin(), demagZ.end(), 0.0);
        std::fill(dipoleX.begin(), dipoleX.end(), 0.0); std::fill(dipoleY.begin(), dipoleY.end(), 0.0); std::fill(dipoleZ.begin(), dipoleZ.end(), 0.0);

        if (systemData->useDipolar)
            Dipolar -> DipolarInteraction1D(mx1, dipoleX);
        if (systemData->useDemagIntense) {
            demagField -> DemagnetisationFieldIntense(demagX, demagY, demagZ, mx1, my1, mz1);
        } else if (systemData->useDemagFft) {
            std::string rkStageName = "2-2";
            demagField -> DemagField1DReal(demagX, demagY, demagZ, mx1, my1, mz1, iteration, rkStageName);
            // if (iteration > 0) {std::cout << "Stage 2" << std::endl; PrintVector(demagZ, false);}
        }

        // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
        for (int site = 1; site <= GV.GetNumSpins(); site++) {

            // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
            int spinLHS = site - 1, spinRHS = site + 1;

            /*double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if (systemData->useDipolar) {
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
            double hxK1 = EffField -> EffectiveFieldX(site, 0, mx1[spinLHS], mx1[site], mx1[spinRHS], dipoleX[site], demagX[site], t0);
            double hyK1 = EffField -> EffectiveFieldY(site, 0, my1[spinLHS], my1[site], my1[spinRHS], dipoleY[site], demagY[site]);
            double hzK1 = EffField -> EffectiveFieldZ(site, 0, mz1[spinLHS], mz1[site], mz1[spinRHS], dipoleZ[site], demagZ[site]);

            // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
            double mxK2 = LLG -> MagneticMomentX(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double myK2 = LLG -> MagneticMomentY(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);
            double mzK2 = LLG -> MagneticMomentZ(site, mx1[site], my1[site], mz1[site], hxK1, hyK1, hzK1);

            mx2[site] = systemData->mx0[site] + systemData->stepsize * mxK2;
            my2[site] = systemData->my0[site] + systemData->stepsize * myK2;
            mz2[site] = systemData->mz0[site] + systemData->stepsize * mzK2;

            if (systemData->shouldTrackMValues) {
                double mIterationNorm = sqrt(pow(mx2[site], 2) + pow(my2[site], 2) + pow(mz2[site], 2));
                if ((systemData->largestMNorm) > (1.0 - mIterationNorm)) { systemData->largestMNorm = (1.0 - mIterationNorm); }
                if (mIterationNorm > 1.00005) {throw std::runtime_error("mag. moments are no longer below <= 1.00005");}
            }
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        systemData->mx0.clear(); systemData->my0.clear(); systemData->mz0.clear();
        mx1.clear(); my1.clear(); mz1.clear();

        childNMData -> SaveDataToFile(mxRK2File, mx2, iteration);
        //childNMData -> SaveDataToFile(myRK2File, my2, iteration);
        //childNMData -> SaveDataToFile(mzRK2File, mz2, iteration);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        systemData->mx0 = mx2; systemData->my0 = my2; systemData->mz0 = mz2;

        if (iteration == systemData->forceStopAtIteration)
            exit(0);

        systemData->totalTime += systemData->stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();
    //myRK2File.close();
    //mzRK2File.close();

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata(true);
    }

    if (systemData->shouldTrackMValues)
        std::cout << "\nMax norm. value of M is: " << systemData->largestMNorm << std::endl;

    // Filename can be copy/pasted from C++ console to Python function's console.
    std::cout << "\n\nFile can be found at:\n\t" << GV.GetFilePath() << GV.GetFileNameBase() << std::endl;
}
void NMMethods::SolveRK2() {
    std::shared_ptr<NMDataHandling> childNMData = std::make_shared<NMDataHandling>(systemData);
    // Uses multiple layers to solve the RK2 midpoint method. See the documentation for more details.

    // Create files to save the data. All files will have (GV.GetFileNameBase()) in them to make them clearly identifiable.
    std::ofstream mxRK2File(GV.GetFilePath() + "rk2_mx_" + GV.GetFileNameBase() + ".csv");
    std::ofstream mxRK2File1(GV.GetFilePath() + "rk2_mx1_" + GV.GetFileNameBase() + ".csv");

    // User information and file header is magnetic-material specific.
    if (systemData->isFm) {
        childNMData -> InformUserOfCodeType("RK2 Midpoint (FM)");
        childNMData -> CreateFileHeader(mxRK2File, "RK2 Midpoint (FM)", false, 0);
        childNMData -> CreateFileHeader(mxRK2File1, "RK2 Midpoint (FM)", false, 1);
    } else if (!systemData->isFm) {
        childNMData -> InformUserOfCodeType("RK2 Midpoint (AFM)");
        childNMData -> CreateFileHeader(mxRK2File, "RK2 Midpoint (AFM)");
        childNMData -> CreateFileHeader(mxRK2File1, "RK2 Midpoint (AFM)");
    }

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata();
    }

    progressbar bar(100);

    std::vector<double> demagX(GV.GetNumSpins() + 2, 0.0), demagY(GV.GetNumSpins() + 2, 0.0), demagZ(GV.GetNumSpins() + 2, 0.0);

    for (int iteration = systemData->iterationStart; iteration <= systemData->iterationEnd; iteration++) {

        if (systemData->iterationEnd >= 100 && iteration % (systemData->iterationEnd / 100) == 0)
            // Doesn't work on Windows due to different compiler. Doesn't work for fewer than 100 iterations
            bar.update();

        _testShockwaveConditions(iteration);

        double t0 = systemData->totalTime;

        for (int layer = 0; layer < systemData->totalLayers; layer++) {
            // RK2 Stage 1. Takes initial conditions as inputs.

            for (int site = 1; site <= systemData->layerTotalSpins[layer]; site++) {
                // Exclude the 0th and last spins as they will always be zero-valued (end, pinned, bound spins)

                // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
                int spinLHS = site - 1, spinRHS = site + 1;

                double mxLHS = systemData->m0Nest[layer][spinLHS][0], mxMID = systemData->m0Nest[layer][site][0], mxRHS = systemData->m0Nest[layer][spinRHS][0];
                double myLHS = systemData->m0Nest[layer][spinLHS][1], myMID = systemData->m0Nest[layer][site][1], myRHS = systemData->m0Nest[layer][spinRHS][1];
                double mzLHS = systemData->m0Nest[layer][spinLHS][2], mzMID = systemData->m0Nest[layer][site][2], mzRHS = systemData->m0Nest[layer][spinRHS][2];

                double dipoleX, dipoleY, dipoleZ;
                if (systemData->useDipolar) {

                    int layer1, layer2;
                    if (layer == 0) {layer1 = 0; layer2 = 1;}
                    else if (layer == 1) {layer1 = 1; layer2 = 0;}

                    if (systemData->debugFunc) {std::cout << "\n\niteration: " << iteration << " | layer: " << layer << " | site: " << site << std::endl;}
                    std::vector<double> dipoleTerms = Dipolar -> DipolarInteractionInterlayer(systemData->m0Nest[layer1], systemData->m0Nest[layer2], site,
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
                double hxK0 = EffField -> EffectiveFieldX(site, layer, mxLHS, mxMID, mxRHS, dipoleX, demagX[site], t0);
                double hyK0 = EffField -> EffectiveFieldY(site, layer, myLHS, myMID, myRHS, dipoleY, demagY[site]);
                double hzK0 = EffField -> EffectiveFieldZ(site, layer, mzLHS, mzMID, mzRHS, dipoleZ, demagZ[site]);

                // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
                double mxK1 = LLG -> MagneticMomentX(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
                double myK1 = LLG -> MagneticMomentY(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);
                double mzK1 = LLG -> MagneticMomentZ(site, layer, mxMID, myMID, mzMID, hxK0, hyK0, hzK0);

                // Find (m0 + k1/2) for each site, which is used in the next stage.
                systemData->m1Nest[layer][site][0] = mxMID + systemData->stepsizeHalf * mxK1;
                systemData->m1Nest[layer][site][1] = myMID + systemData->stepsizeHalf * myK1;
                systemData->m1Nest[layer][site][2] = mzMID + systemData->stepsizeHalf * mzK1;
            }
        }

        for (int layer = 0; layer < systemData->totalLayers; layer++) {
            // RK2 Stage 2. Takes (m0 + k1/2) as inputs.
            for (int site = 1; site <= systemData->layerTotalSpins[layer]; site++) {

                // Relative to the current site (site); site to the left (LHS); site to the right (RHS)
                int spinLHS = site - 1, spinRHS = site + 1;

                double mxLHS = systemData->m1Nest[layer][spinLHS][0], mxMID = systemData->m1Nest[layer][site][0], mxRHS = systemData->m1Nest[layer][spinRHS][0];
                double myLHS = systemData->m1Nest[layer][spinLHS][1], myMID = systemData->m1Nest[layer][site][1], myRHS = systemData->m1Nest[layer][spinRHS][1];
                double mzLHS = systemData->m1Nest[layer][spinLHS][2], mzMID = systemData->m1Nest[layer][site][2], mzRHS = systemData->m1Nest[layer][spinRHS][2];

                double dipoleX, dipoleY, dipoleZ;
                if (systemData->useDipolar) {

                    int layer1, layer2;
                    if (layer == 0) {layer1 = 0; layer2 = 1;}
                    else if (layer == 1) {layer1 = 1; layer2 = 0;}

                    int debugCounter = 0;  // To make sure debug outputs only occur during the first RK2 stage, not this second stage
                    if (systemData->debugFunc) { systemData->debugFunc = false; debugCounter++; }
                    std::vector<double> dipoleTerms = Dipolar -> DipolarInteractionInterlayer(systemData->m1Nest[layer1], systemData->m1Nest[layer2], site,
                                                                                   layer1, layer2);
                    if (debugCounter > 0) { systemData->debugFunc = true; }

                    dipoleX = dipoleTerms[0];
                    dipoleY = dipoleTerms[1];
                    dipoleZ = dipoleTerms[2];
                } else {
                    dipoleX = 0;
                    dipoleY = 0;
                    dipoleZ = 0;
                }
                // Calculations for the effective field (H_eff), coded as symbol 'h', components of the target site
                double hxK1 = EffField -> EffectiveFieldX(site, layer, mxLHS, mxMID, mxRHS, dipoleX, demagX[site], t0);
                double hyK1 = EffField -> EffectiveFieldY(site, layer, myLHS, myMID, myRHS, dipoleY, demagY[site]);
                double hzK1 = EffField -> EffectiveFieldZ(site, layer, mzLHS, mzMID, mzRHS, dipoleZ, demagZ[site]);

                // RK2 K-value calculations for the magnetic moment, coded as symbol 'm', components of the target site
                double mxK2 = LLG -> MagneticMomentX(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
                double myK2 = LLG -> MagneticMomentY(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);
                double mzK2 = LLG -> MagneticMomentZ(site, layer, mxMID, myMID, mzMID, hxK1, hyK1, hzK1);

                systemData->m2Nest[layer][site][0] = systemData->m0Nest[layer][site][0] + systemData->stepsize * mxK2;
                systemData->m2Nest[layer][site][1] = systemData->m0Nest[layer][site][1] + systemData->stepsize * myK2;
                systemData->m2Nest[layer][site][2] = systemData->m0Nest[layer][site][2] + systemData->stepsize * mzK2;

                if (systemData->shouldTrackMValues) {
                    double mIterationNorm = sqrt(
                            pow(systemData->m2Nest[layer][site][0], 2) + pow(systemData->m2Nest[layer][site][1], 2) + pow(systemData->m2Nest[layer][site][2], 2));
                    if ((systemData->largestMNormMulti[layer]) > (1.0 - mIterationNorm)) { systemData->largestMNormMulti[layer] = (1.0 - mIterationNorm); }
                }
            }
        }

        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */

        childNMData -> SaveDataToFileMultilayer(mxRK2File, systemData->m2Nest[0], iteration, 0);
        childNMData -> SaveDataToFileMultilayer(mxRK2File1, systemData->m2Nest[1], iteration, 1);

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        systemData->m0Nest = systemData->m2Nest;

        if (iteration == systemData->forceStopAtIteration)
            exit(0);

        systemData->totalTime += systemData->stepsize;
    }// Final line of RK2 solver for all iterations. Everything below here occurs after RK2 method is complete

    // Ensures files are closed; sometimes are left open if the writing process above fails
    mxRK2File.close();

    if (GV.GetEmailWhenCompleted()) {
        childNMData -> CreateMetadata(true);
    }

    if (systemData->shouldTrackMValues) {
        std::cout << "\nMax norm. values of M are: ";
        for (int i = 0; i < systemData->largestMNormMulti.size(); i++) {
            std::cout << "Layer " << i << ": " << systemData->largestMNormMulti[i] << " | ";
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