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

SolversSuperClass::SolversSuperClass(std::shared_ptr<SimulationParameters> sharedSimParams,
                                     std::shared_ptr<SimulationStates> sharedSimStates, 
                                     std::shared_ptr<SimulationFlags> sharedSimFlags)  : 
                                     
    simParams(std::move(sharedSimParams)), simStates(std::move(sharedSimStates)), simFlags(std::move(sharedSimFlags)) {}

SolversSuperClass::~SolversSuperClass() {
    // Memory automatically deallocated by unique_ptr
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createSimulationInstance(std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                               std::shared_ptr<SimulationStates> sharedSimStates, 
                                                                               std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<NMInitialisation>(sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createConfigurationInstance(std::shared_ptr<SimulationParameters> sharedSimParams, 
                                                                                  std::shared_ptr<SimulationStates> sharedSimStates, 
                                                                                  std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<NMConfiguration>(sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createDataHandlingInstance(std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                                 std::shared_ptr<SimulationStates> sharedSimStates,
                                                                                 std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<NMDataHandling>(sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createMethodsInstance(std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                            std::shared_ptr<SimulationStates> sharedSimStates,
                                                                            std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<NMMethods>(sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SimulationParameters> SolversSuperClass::getSimulationParamsContainer() {
    // Lazy initialisation. Left more as a reminder of how to do it than anything else
    static std::shared_ptr<SimulationParameters> sharedSimParams = std::make_shared<SimulationParameters>();
    return sharedSimParams;
}

std::shared_ptr<SimulationStates> SolversSuperClass::getSimulationStatesContainer() {
    // Lazy initialisation. Left more as a reminder of how to do it than anything else
    static std::shared_ptr<SimulationStates> sharedSimStates = std::make_shared<SimulationStates>();
    return sharedSimStates;
}

std::shared_ptr<SimulationFlags> SolversSuperClass::getSimulationFlagsContainer() {
    // Lazy initialisation. Left more as a reminder of how to do it than anything else
    static std::shared_ptr<SimulationFlags> sharedSimFlags = std::make_shared<SimulationFlags>();
    return sharedSimFlags;
}

/*
std::pair<double, double> SolversSuperClass::testInstance() {
    auto instance = std::dynamic_pointer_cast<NMInitialisation>(createSimulationInstance(sharedSimParams));
    double before = instance->getsharedSimParamsContainer()->ambientTemperature; // Get value of var1 before modification

    instance->testModifyingDouble(100); // Modify var1 using ChildClass method

    std::cout << "This works because simParams is protected, and the parent can still access it:" << instance->simParams->iterEndShock << std::endl;

    double after = instance->getsharedSimParamsContainer()->ambientTemperature; // Get value of var1 after modification
    return {before, after};
}
 */