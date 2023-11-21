//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMMETHODS_H
#define SPINCHAINS_NMMETHODS_H

// C++ User Libraries (Parent Class)
#include <utility>
#include "../libs/progressbar.hpp"
#include "SolversSuperClass.h"

// C++ User Libraries (Sibling Classes)
#include "NMDataHandling.h"

// C++ User Libraries (Interface)
#include "INMMethods.h"

// C++ User Libraries (Components)
#include "DemagField.h"
#include "EffectiveField.h"
#include "LLG.h"
#include "DipolarField.h"

class NMMethods :  public SolversSuperClass, public INMMethods {
private:
    DemagnetisationFields demagField;
    EffectiveField effectiveField;
    MagnetisationDynamics llg;
    DipolarInteractions dipolarField;
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
    NMMethods(std::shared_ptr<SimulationParameters> paramsData,
              std::shared_ptr<SimulationStates> sharedSimStates,
              std::shared_ptr<SimulationFlags> sharedSimFlags);
    ~NMMethods() = default;
public:
    void performInitialisation() override { runMethod(); };

    void                runMethod() override;
};


#endif //SPINCHAINS_NMMETHODS_H
