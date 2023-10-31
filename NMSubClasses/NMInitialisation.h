//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMINITIALISATION_H
#define SPINCHAINS_NMINITIALISATION_H

#include "../NMSuperClassTest.h"

class NMInitialisation: public NMSuperClassTest {
private:
    // ####################################            Define Private Variables            ###################################

    double              _recordingInterval;
    int                 _layerOfInterest;
    const double        _BOHR_MAGNETON = 9.274e-24;                               // Bohr magneton in Am^{2} (equiv. to J T^{-1})




    // ####################################            Define Private Functions            ###################################

    // Description missing
    void                _setSimulationFlags();

    // Description missing
    void                _setSimulationParameters();

    // Description missing
    void                _generateRemainingParameters();

    // Description missing
    void                _setMaterialParameters();

    // Description missing
    void                _guardClauses();

protected:
    // ####################################            Define Protected Variables            ###################################

    // ####################################            Define Protected Functions            ###################################


public:
    // ####################################            Define Public Variables            ###################################

    // ####################################            Define Public Functions            ###################################

    // Initialisation constructor
    void                Initialise();
};


#endif //SPINCHAINS_NMINITIALISATION_H
