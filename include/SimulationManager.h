//
// Created by Cameron McEleney on 02/02/2024.
//

#ifndef SPINCHAINS_SIMULATIONMANAGER_H
#define SPINCHAINS_SIMULATIONMANAGER_H

#include <string>

class SimulationManager {
public:
    std::string             outputFileID;
    bool                    hasNumericSuffix = true;
    bool                    resetSimState = false;
    bool                    massProduce = false;
    bool                    hasFirstRunOccurred = false;
};

#endif //SPINCHAINS_SIMULATIONMANAGER_H
