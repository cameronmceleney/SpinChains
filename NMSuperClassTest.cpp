//
// Created by Cameron McEleney on 31/10/2023.
//

// C++ User Libraries (Parent)
#include "NMSuperClassTest.h"

// C++ User Libraries (Children)
#include "NMSubClasses/NMInitialisation.h"

NMSuperClassTest::NMSuperClassTest(std::shared_ptr<SharedVariableHolder> data) : sharedData(data) {}

NMSuperClassTest::NMSuperClassTest() {
    //
}

NMSuperClassTest::~NMSuperClassTest() {
    // Memory automatically deallocated by unique_ptr
}

void NMSuperClassTest::callInitialise() {
    // Default implementation or no implementation
}

void NMSuperClassTest::executeDerivedMethod() {
    std::cout << sharedVariables.ambientTemperature << std::endl;
    std::shared_ptr<NMInitialisation> derivedPtr = std::make_shared<NMInitialisation>(sharedData);
    derivedPtr -> callInitialise();
    std::cout << sharedVariables.ambientTemperature << std::endl;
}