//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSSUPERCLASS_H
#define SPINCHAINS_SOLVERSSUPERCLASS_H

// C++ Standard Library
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <memory>

// C++ User Libraries (General)
#include "SystemDataContainer.h"
#include "../libs/linspace.h"

class SolversSuperClass {
private:
    //SystemDataContainer sharedVariables;#

protected:
    // Getter for child classes to access sharedVariables
    std::shared_ptr<SystemDataContainer> simState;

public:
    //SolversSuperClass(std::shared_ptr<SystemDataContainer> simState);
    SolversSuperClass(std::shared_ptr<SystemDataContainer> data) : simState(std::move(data)) {}
    virtual ~SolversSuperClass(); // Destructor (though not strictly necessary with smart pointers)

public:
    virtual void performInitialisation() = 0;
    static std::shared_ptr<SolversSuperClass> createSimulationInstance(std::shared_ptr<SystemDataContainer> sharedData);
    static std::shared_ptr<SolversSuperClass> createConfigurationInstance(std::shared_ptr<SystemDataContainer> sharedData);
    static std::shared_ptr<SolversSuperClass> createDataHandlingInstance(std::shared_ptr<SystemDataContainer> sharedData);
    static std::shared_ptr<SolversSuperClass> createMethodsInstance(std::shared_ptr<SystemDataContainer> sharedData);

    static std::shared_ptr<SystemDataContainer> getSharedDataContainer();

    static std::pair<double, double> testInstance();
};


#endif //SPINCHAINS_SOLVERSSUPERCLASS_H
