//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSIMPLEMENTATION_H
#define SPINCHAINS_SOLVERSIMPLEMENTATION_H

// C++ Corresponding Interface
#include "InterfaceSolversImplementation.h"

// C++ Standard Libraries
#include <string>

// C++ User Libraries (General)
#include "GlobalVariables.h"
#include "../libs/progressbar.hpp"

// C++ User Libraries (Parent Class)
#include "SolversSuperClass.h"
// C++ User Libraries (Sibling Classes)
#include "SolversDataHandling.h"

// C++ User Libraries (Class' Components)
#include "DemagnetisationFields.h"
#include "EffectiveFields.h"
#include "MagnetisationDynamics.h"
#include "DipolarFields.h"

class SolversImplementation :  public SolversSuperClass, public iSolversImplementation {
private:
    DemagnetisationFields demagField;
    EffectiveFields effectiveField;
    MagnetisationDynamics llg;
    DipolarFields dipolarField;
private:
    // Description missing
    void                _testShockwaveConditions(double iteration) override;
            // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2Classic() override;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2() override;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    // void                SolveRK2Test();

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
