//
// Created by Cameron McEleney on 31/10/2023.
//
#pragma once

#ifndef SPINCHAINS_EFFECTIVEFIELDS_H
#define SPINCHAINS_EFFECTIVEFIELDS_H

// C++ Standard Library
#include <vector>

// C++ User Libraries (General)
#include "CommonLibs.h"

// C++ User Libraries (General)
#include "GlobalVariables.h"

// C++ User Libraries (Containers)
#include "SimulationParameters.h"
#include "SimulationStates.h"
#include "SimulationFlags.h"

class EffectiveFields {
private:
    SimulationParameters *_simParams; // Non-owning pointer to SimulationParameters
    SimulationStates *_simStates;
    SimulationFlags *_simFlags;
public:
    explicit EffectiveFields( SimulationParameters *sharedSimParams,
                              SimulationStates *sharedSimStates,
                              SimulationFlags *sharedSimFlags );

    ~EffectiveFields() = default;

public:
    double EffectiveFieldXClassic( const int &site, const int &layer, const double &mxLHS, const double &mxMID,
                                   const double &mxRHS, const double &dipoleTerm, const double &demagTerm,
                                   const double &dmiTerm, const double &current_time );

    // Description missing
    double EffectiveFieldYClassic( const int &site, const int &layer, const double &myLHS, const double &myMID,
                                   const double &myRHS, const double &dipoleTerm, const double &demagTerm,
                                   const double &dmiTerm );

    // Description missing
    double EffectiveFieldZClassic( const int &site, const int &layer, const double &mzLHS, const double &mzMID,
                                   const double &mzRHS, const double &dipoleTerm,
                                   const double &demagTerm );

    void EffectiveFieldXTest( const int &currentSite, const int &layer, const std::vector<double> &mxTermsIn,
                              std::vector<double> &effectiveFieldXTermsOut, const double &current_time );

    // Description missing
    void EffectiveFieldYTest( const int &currentSite, const int &layer, const std::vector<double> &myTermsIn,
                              std::vector<double> &effectiveFieldYTermsOut );

    // Description missing
    void EffectiveFieldZTest( const int &currentSite, const int &layer, const std::vector<double> &mzTermsIn,
                              std::vector<double> &effectiveFieldZTermsOut );

    void EffectiveFieldsCombinedTest( const int &currentSite, const int &layer, const std::vector<double> &mxTermsIn,
                                      const std::vector<double> &myTermsIn, const std::vector<double> &mzTermsIn,
                                      std::vector<double> &effectiveFieldXTermsOut,
                                      std::vector<double> &effectiveFieldYTermsOut,
                                      std::vector<double> &effectiveFieldZTermsOut, const double &currentTime );

    std::array<double, 3>
    EffectiveFieldsCombinedTestExOnly( const int &currentSite, const int &layer,
                                       const std::vector<double> &mxTermsIn,
                                       const std::vector<double> &myTermsIn,
                                       const std::vector<double> &mzTermsIn );

    std::array<double, 3> EffectiveFieldsCombinedTestDriveOnly( const int &currentSite, const int &layer,
                                                                const std::vector<double> &mxTermsIn,
                                                                const std::vector<double> &myTermsIn,
                                                                const std::vector<double> &mzTermsIn,
                                                                const double &currentTime );

    bool isSiteDriven( const int &site );

    struct exchangeOut {
        double x;
        double y;
        double z;

        exchangeOut() : x(0.0), y(0.0), z(0.0) {}

        exchangeOut( double updateX, double updateY, double updateZ ) : x(updateX), y(updateY), z(updateZ) {}
    };
};


#endif //SPINCHAINS_EFFECTIVEFIELDS_H
