//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMINITIALISATION_H
#define SPINCHAINS_NMINITIALISATION_H

#include "SolversSuperClass.h"

class NMInitialisation: public SolversSuperClass {
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

    void                Initialise();


protected:
    //

public:
    NMInitialisation(std::shared_ptr<SystemDataContainer> data) : SolversSuperClass(data){};
    void                testModifyingDouble(double  newValue);
    void                performInitialisation() override { Initialise(); };
};


#endif //SPINCHAINS_NMINITIALISATION_H
