//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMMETHODS_H
#define SPINCHAINS_NMMETHODS_H

// C++ User Libraries (Parent Class)
#include <utility>

#include "../include/NMSuperClassTest.h"
// C++ User Libraries (Sibling Classes)
#include "NMDataHandling.h"
#include "INMMethods.h"
// C++ User Libraries (Child Classes)
#include "IDemagField.h"
#include "IDipolarField.h"
#include "IEffectiveField.h"
#include "ILLG.h"

class NMMethods :  public NMSuperClassTest, public INMMethods {
private:
    // ####################################            Private Instances            ###################################
    //NMDataHandling* NMData;

    std::shared_ptr<IDemagField> demagField;
    std::shared_ptr<IEffectiveField> EffField;
    std::shared_ptr<ILLG> LLG;
    std::shared_ptr<IDipolarField> Dipolar;
    // Description missing
    void                _testShockwaveConditions(double iteration) override;
            // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2Classic() override;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2() override;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    // void                SolveRK2Test();

public:
    NMMethods(std::shared_ptr<SystemDataContainer> data);
    void performInitialisation() override { runMethod(); };

public:
    void                runMethod() override;
};


#endif //SPINCHAINS_NMMETHODS_H
