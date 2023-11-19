//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMMETHODS_H
#define SPINCHAINS_NMMETHODS_H

// C++ User Libraries (Parent Class)
#include "../include/NMSuperClassTest.h"

// C++ User Libraries (Sibling Classes)
#include "NMDataHandling.h"

// C++ User Libraries (Child Classes)
class DemagnetisationFields;
class EffectiveField;
class MagnetisationDynamics;
class DebuggingTools;
class DipolarInteractions;

class NMMethods:  public NMSuperClassTest {
private:
    // ####################################            Private Instances            ###################################
    //NMDataHandling* NMData;

    DemagnetisationFields* Demag;
    EffectiveField* EffField;
    MagnetisationDynamics* LLG;
    DebuggingTools* Debug;
    DipolarInteractions* Dipolar;
    // ####################################            Define Private Variables            ###################################

    // ####################################            Define Private Methods            ###################################


    // Description missing
    void                _testShockwaveConditions(double iteration);
            // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2Classic();

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2();

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    // void                SolveRK2Test();



protected:
    // ####################################            Protected Instances            ###################################

    // ####################################            Define Protected Variables            ###################################

    // ####################################            Define Protected Functions            ###################################


public:
    // ####################################            Define Public Variables            ###################################

    // ####################################            Define Public Methods            ###################################
    NMMethods(std::shared_ptr<SystemDataContainer> data) : NMSuperClassTest(data){};
    void performInitialisation() override { runMethod(); };

public:
    void                runMethod();
};


#endif //SPINCHAINS_NMMETHODS_H
