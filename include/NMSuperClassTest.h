//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMSUPERCLASSTEST_H
#define SPINCHAINS_NMSUPERCLASSTEST_H

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
#include "../NMSubClasses/SystemDataContainer.h"
#include "linspace.h"

class NMSuperClassTest {
private:
    //SystemDataContainer sharedVariables;#

protected:
    // Getter for child classes to access sharedVariables
    std::shared_ptr<SystemDataContainer> simState;

public:
    //NMSuperClassTest(std::shared_ptr<SystemDataContainer> simState);
    NMSuperClassTest(std::shared_ptr<SystemDataContainer> data) : simState(std::move(data)) {}
    virtual ~NMSuperClassTest(); // Destructor (though not strictly necessary with smart pointers)

public:
    virtual void performInitialisation() = 0;
    static std::shared_ptr<NMSuperClassTest> createSimulationInstance(std::shared_ptr<SystemDataContainer> sharedData);
    static std::shared_ptr<NMSuperClassTest> createConfigurationInstance(std::shared_ptr<SystemDataContainer> sharedData);
    static std::shared_ptr<NMSuperClassTest> createDataHandlingInstance(std::shared_ptr<SystemDataContainer> sharedData);
    static std::shared_ptr<NMSuperClassTest> createMethodsInstance(std::shared_ptr<SystemDataContainer> sharedData);

    static std::shared_ptr<SystemDataContainer> getSharedDataContainer();

    static std::pair<double, double> testInstance();
};


#endif //SPINCHAINS_NMSUPERCLASSTEST_H
