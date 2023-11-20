//
// Created by Cameron McEleney on 31/10/2023.
//

#include <utility>

// C++ User Libraries (Parent)
#include "../include/NMSuperClassTest.h"

// C++ User Libraries (Children)
#include "../NMSubClasses/NMInitialisation.h"
#include "../NMSubClasses/NMConfiguration.h"
#include "../NMSubClasses/NMDataHandling.h"
#include "../NMSubClasses/NMMethods.h"

NMSuperClassTest::~NMSuperClassTest() {
    // Memory automatically deallocated by unique_ptr
}

std::shared_ptr<NMSuperClassTest> NMSuperClassTest::createSimulationInstance() {
    // Return instance of the child
    auto sharedData = getSharedDataContainer();
    return std::make_shared<NMInitialisation>(sharedData);
}

std::shared_ptr<NMSuperClassTest> NMSuperClassTest::createConfigurationInstance() {
    // Return instance of the child
    auto sharedData = getSharedDataContainer();
    return std::make_shared<NMConfiguration>(sharedData);
}

std::shared_ptr<NMSuperClassTest> NMSuperClassTest::createDataHandlingInstance() {
    // Return instance of the child
    auto sharedData = getSharedDataContainer();
    return std::make_shared<NMDataHandling>(sharedData);
}

std::shared_ptr<NMSuperClassTest> NMSuperClassTest::createMethodsInstance() {
    // Return instance of the child
    auto sharedData = getSharedDataContainer();
    return std::make_shared<NMMethods>(sharedData);
}

std::shared_ptr<SystemDataContainer> NMSuperClassTest::getSharedDataContainer() {
    // Lazy initialisation. Left more as a reminder of how to do it than anything else
    static std::shared_ptr<SystemDataContainer> sharedData = std::make_shared<SystemDataContainer>();
    return sharedData;
}

std::pair<double, double> NMSuperClassTest::testInstance() {
    auto instance = std::dynamic_pointer_cast<NMInitialisation>(createSimulationInstance());
    double before = instance->getSharedDataContainer()->ambientTemperature; // Get value of var1 before modification

    instance->testModifyingDouble(100); // Modify var1 using ChildClass method

    std::cout << "This works because systemData is protected, and the parent can still access it:" << instance->systemData->iterEndShock << std::endl;

    double after = instance->getSharedDataContainer()->ambientTemperature; // Get value of var1 after modification
    return {before, after};
}