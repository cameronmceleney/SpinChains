//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSIMPLEMENTATION_H
#define SPINCHAINS_SOLVERSIMPLEMENTATION_H

// C++ Corresponding Interface
#include "InterfaceSolversImplementation.h"

// C++ Standard Libraries
#include <string>

// C++ Third Party Library
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

// C++ User Libraries (General)
#include "GlobalVariables.h"
#include "../libs/progressbar.hpp"

// C++ User Libraries (Parent Class)
#include "SolversSuperClass.h"
// C++ User Libraries (Sibling Classes)
#include "SolversDataHandling.h"

// C++ User Libraries (Class' Components)
#include "DemagnetisationFields.h"
#include "DzyaloshinskiiMoriyaInteraction.h"
#include "EffectiveFields.h"
#include "MagnetisationDynamics.h"
#include "DipolarFields.h"

class SolversImplementation :  public SolversSuperClass, public iSolversImplementation {
private:
    DemagnetisationFields demagField;
    DzyaloshinskiiMoriyaInteraction dmInteraction;
    EffectiveFields effectiveField;
    MagnetisationDynamics llg;
    DipolarFields dipolarField;
private:
    /**
     * Description missing
     */
    void                _testShockwaveConditions(double iteration) override;
    /**
     * Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method. Original
     * implementation of the RK2 method; used for testing purposes
     */
    void                SolveRK2Classic() override;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2() override;

    /**
     * Evaluate the given system using the RK2 method which has been configured to enabled
     * parallelisation through Intel's oneAPI TBB
     */
    void                RK2Parallel() override;

    /**
     * Description missing
     */
    void RK2StageMultithreaded(const std::vector<double> &mxIn, const std::vector<double> &myIn,
                               const std::vector<double> &mzIn, std::vector<double> &mxOut,
                               std::vector<double> &myOut, std::vector<double> &mzOut,
                               std::vector<double> &demagX, std::vector<double> &demagY,
                               std::vector<double> &demagZ, std::vector<double> &dipoleX,
                               std::vector<double> &dipoleY, std::vector<double> &dipoleZ,
                               std::vector<double> &dmiX, std::vector<double> &dmiY,
                               std::vector<double> &dmiZ, double &currentTime, double &stepsize,
                               int &iteration, std::string rkStage);
public:
    SolversImplementation(std::shared_ptr<SimulationParameters> paramsData,
              std::shared_ptr<SimulationStates> sharedSimStates,
              std::shared_ptr<SimulationFlags> sharedSimFlags);
    ~SolversImplementation() = default;
public:
    void performInitialisation() override { runMethod(); };

    void                runMethod() override;
};


#endif //SPINCHAINS_SOLVERSIMPLEMENTATION_H
