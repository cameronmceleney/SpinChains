//
// Created by Cameron McEleney on 31/10/2023.
//

#include <utility>

// C++ User Libraries (Parent)
#include "../include/NMSuperClassTest.h"

// C++ User Libraries (Children)
#include "../NMSubClasses/NMInitialisation.h"


NMSuperClassTest::~NMSuperClassTest() {
    // Memory automatically deallocated by unique_ptr
}

std::shared_ptr<NMSuperClassTest> NMSuperClassTest::createSimulationInstance() {
    // Return instance of the child
    return std::make_shared<NMInitialisation>();
}

std::shared_ptr<SystemDataContainer> NMSuperClassTest::getSystemData() {
    // Lazy initialisation. Left more as a reminder of how to do it than anything else
    if (!systemData)
        systemData = std::make_shared<SystemDataContainer>();

    return systemData;
}

std::pair<double, double> NMSuperClassTest::testInstance() {
    auto instance = std::dynamic_pointer_cast<NMInitialisation>(createSimulationInstance());
    double before = instance->getSystemData()->ambientTemperature; // Get value of var1 before modification

    instance->testModifyingDouble(100); // Modify var1 using ChildClass method

    std::cout << "This works because systemData is protected, and the parent can still access it:" << instance->systemData->iterEndShock << std::endl;

    double after = instance->getSystemData()->ambientTemperature; // Get value of var1 after modification
    return {before, after};
}