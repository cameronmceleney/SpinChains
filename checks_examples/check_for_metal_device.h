//
// check_for_metal_device.h
// Created by Cameron Aidan McEleney on 08/03/2024.
//

#ifndef SPINCHAINS_CHECK_FOR_METAL_DEVICE_H
#define SPINCHAINS_CHECK_FOR_METAL_DEVICE_H

#include "CommonTools.h"

#include <memory>
#include <string>
#include <vector>

// Forward declarations for types e.g.
//namespace KeyFlags {
//    class FullDebug;
//}

// Declaration of CoreParams structs as before
class CoffeeExample: public CommonTools {
public:

    CoffeeExample();
    ~CoffeeExample(); // Implement to define the destructor where Impl is a complete type

    // Public interface as before
    void vectorAddition( const int &numElements, const std::string &kernelFunctionName, const bool &fullDebug );
    // More public methods...

private:
    class CoffeeExampleImpl; // Forward declaration
    std::unique_ptr<CoffeeExampleImpl> pCoffeeExampleImpl; // Opaque pointer


};

#endif //SPINCHAINS_CHECK_FOR_METAL_DEVICE_H
