//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMINITIALISATION_H
#define SPINCHAINS_NMINITIALISATION_H

#include "../include/NMSuperClassTest.h"

class NMInitialisation: public NMSuperClassTest {
private:
    // ####################################            Define Private Variables            ###################################

    double              _recordingInterval;
    int                 _layerOfInterest;
    const double        _BOHR_MAGNETON = 9.274e-24;


    // ####################################            Define Private Functions            ###################################

    void                _setSimulationFlags();
    void                _setSimulationParameters();
    void                _generateRemainingParameters();
    void                _setMaterialParameters();
    void                _guardClauses();

protected:
    //

public:
    //NMInitialisation(std::shared_ptr<SystemDataContainer> data);
    void                initialiseSimulation() override;
    void                testModifyingDouble(double  newValue);
    void                Initialise();
};


#endif //SPINCHAINS_NMINITIALISATION_H
