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

// C++ Third Party Libraries
extern "C" {
    #include <fftw3.h>
}

// C++ User Libraries (General)
#include "linspace.h"
#include "SpinChainEigenSolverClass.h"
#include "../src/CommonLibs.h"
#include "progressbar.hpp"
#include "../NMSubClasses/SystemDataContainer.h"

class NMSuperClassTest {
private:
    //SystemDataContainer sharedVariables;#

protected:
    // Getter for child classes to access sharedVariables
    std::shared_ptr<SystemDataContainer> systemData;

public:
    //NMSuperClassTest(std::shared_ptr<SystemDataContainer> systemData);
    explicit NMSuperClassTest(std::shared_ptr<SystemDataContainer> data) : systemData(std::move(data)) {}
    virtual ~NMSuperClassTest(); // Destructor (though not strictly necessary with smart pointers)

public:
    virtual void performInitialisation() = 0;
    static std::shared_ptr<NMSuperClassTest> createSimulationInstance();
    static std::shared_ptr<NMSuperClassTest> createConfigurationInstance();
    static std::shared_ptr<NMSuperClassTest> createDataHandlingInstance();
    static std::shared_ptr<NMSuperClassTest> createMethodsInstance();

    static std::shared_ptr<SystemDataContainer> getSharedDataContainer();

    static std::pair<double, double> testInstance();
};


#endif //SPINCHAINS_NMSUPERCLASSTEST_H
