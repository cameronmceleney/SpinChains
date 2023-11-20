//
// Created by Cameron McEleney on 31/10/2023.
//

#include <utility>

// C++ User Libraries (Parent)
#include "../include/SolversSuperClass.h"

// C++ User Libraries (Children)
#include "../include/NMInitialisation.h"
#include "../include/NMConfiguration.h"
#include "../include/NMDataHandling.h"
#include "../include/NMMethods.h"

SolversSuperClass::~SolversSuperClass() {
    // Memory automatically deallocated by unique_ptr
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createSimulationInstance(std::shared_ptr<SystemDataContainer> sharedData) {
    // Return instance of the child
    //auto sharedData = getSharedDataContainer();
    return std::make_shared<NMInitialisation>(sharedData);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createConfigurationInstance(std::shared_ptr<SystemDataContainer> sharedData) {
    // Return instance of the child
    //auto sharedData = getSharedDataContainer();
    return std::make_shared<NMConfiguration>(sharedData);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createDataHandlingInstance(std::shared_ptr<SystemDataContainer> sharedData) {
    // Return instance of the child
    //auto sharedData = getSharedDataContainer();
    return std::make_shared<NMDataHandling>(sharedData);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createMethodsInstance(std::shared_ptr<SystemDataContainer> sharedData) {
    // Return instance of the child
    //auto sharedData = getSharedDataContainer();
    return std::make_shared<NMMethods>(sharedData);
}

std::shared_ptr<SystemDataContainer> SolversSuperClass::getSharedDataContainer() {
    // Lazy initialisation. Left more as a reminder of how to do it than anything else
    static std::shared_ptr<SystemDataContainer> sharedData = std::make_shared<SystemDataContainer>();
    return sharedData;
}

/*
std::pair<double, double> SolversSuperClass::testInstance() {
    auto instance = std::dynamic_pointer_cast<NMInitialisation>(createSimulationInstance(sharedData));
    double before = instance->getSharedDataContainer()->ambientTemperature; // Get value of var1 before modification

    instance->testModifyingDouble(100); // Modify var1 using ChildClass method

    std::cout << "This works because simState is protected, and the parent can still access it:" << instance->simState->iterEndShock << std::endl;

    double after = instance->getSharedDataContainer()->ambientTemperature; // Get value of var1 after modification
    return {before, after};
}
 */