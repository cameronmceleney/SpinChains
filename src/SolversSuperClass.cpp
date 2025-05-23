//
// Created by Cameron McEleney on 31/10/2023.
//

// Corresponding header
#include "../include/SolversSuperClass.h"

// C++ User Libraries (Class' Children)
#include"../include/SolversManager.h"
#include "../include/SolversInitialisation.h"
#include "../include/SolversConfiguration.h"
#include "../include/SolversDataHandling.h"
#include "../include/SolversImplementation.h"  // not used due to inclusion in SolversDataHandling (for now; that'll change in the future)

SolversSuperClass::SolversSuperClass(std::shared_ptr<SimulationManager> sharedSimManager,
                                     std::shared_ptr<SimulationParameters> sharedSimParams,
                                     std::shared_ptr<SimulationStates> sharedSimStates, 
                                     std::shared_ptr<SimulationFlags> sharedSimFlags)  : 
                                     
    simManager(std::move(sharedSimManager)), simParams(std::move(sharedSimParams)), simStates(std::move(sharedSimStates)), simFlags(std::move(sharedSimFlags)) {}

SolversSuperClass::~SolversSuperClass() {
    // Memory automatically deallocated by unique_ptr
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createSimulationManager(std::shared_ptr<SimulationManager> sharedSimMananger,
                                                                              std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                              std::shared_ptr<SimulationStates> sharedSimStates,
                                                                              std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<SolversManager>(sharedSimMananger, sharedSimParams, sharedSimStates, sharedSimFlags);
}


std::shared_ptr<SolversSuperClass> SolversSuperClass::createSimulationInstance(std::shared_ptr<SimulationManager> sharedSimMananger,
                                                                               std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                               std::shared_ptr<SimulationStates> sharedSimStates, 
                                                                               std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<SolversInitialisation>(sharedSimMananger, sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createConfigurationInstance(std::shared_ptr<SimulationManager> sharedSimMananger,
                                                                                  std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                                  std::shared_ptr<SimulationStates> sharedSimStates, 
                                                                                  std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<SolversConfiguration>(sharedSimMananger, sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createDataHandlingInstance(std::shared_ptr<SimulationManager> sharedSimMananger,
                                                                                 std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                                 std::shared_ptr<SimulationStates> sharedSimStates,
                                                                                 std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<SolversDataHandling>(sharedSimMananger, sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SolversSuperClass> SolversSuperClass::createMethodsInstance(std::shared_ptr<SimulationManager> sharedSimMananger,
                                                                            std::shared_ptr<SimulationParameters> sharedSimParams,
                                                                            std::shared_ptr<SimulationStates> sharedSimStates,
                                                                            std::shared_ptr<SimulationFlags> sharedSimFlags) {
    // Return instance of the child
    //auto sharedSimParams = getsharedSimParamsContainer();
    return std::make_shared<SolversImplementation>(sharedSimMananger, sharedSimParams, sharedSimStates, sharedSimFlags);
}

std::shared_ptr<SimulationManager> SolversSuperClass::getSimulationManagerContainer() {
    // Lazy initialisation. Left more as a reminder of how to do it than anything else
    static std::shared_ptr<SimulationManager> sharedSimManager = std::make_shared<SimulationManager>();
    return sharedSimManager;
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
    auto instance = std::dynamic_pointer_cast<SolversInitialisation>(createSimulationInstance(sharedSimParams));
    double before = instance->getsharedSimParamsContainer()->ambientTemperature; // Get value of var1 before modification

    instance->testModifyingDouble(100); // Modify var1 using ChildClass method

    std::cout << "This works because simParams is protected, and the parent can still access it:" << instance->simParams->risingTimeEndAtIteration << std::endl;

    double after = instance->getsharedSimParamsContainer()->ambientTemperature; // Get value of var1 after modification
    return {before, after};
}
 */