//
// Created by Cameron Aidan McEleney on 03/04/2024.
//

#ifndef SPINCHAINS_LOCALMETALRESOURCES_H
#define SPINCHAINS_LOCALMETALRESOURCES_H

#include "../include/SimulationFlags.h"
#include "../include/SimulationParameters.h"
#include "../include/SimulationStates.h"

#include <iostream>
#include <array>
#include <Metal/Metal.hpp>
struct MetalResources {
    // Manual autopool of resources
    MTL::Device *device = nullptr;
    MTL::CommandQueue *commandQueuePrimary = nullptr;
    MTL::ComputePipelineState *computePrimaryPipelineState = nullptr;
    NS::Error *error = nullptr;
    NS::String *libraryPath = nullptr;
    MTL::Library *library = nullptr;
    MTL::Function *kernelFunction = nullptr;
    MTL::CommandBuffer *computeCommandBuffer = nullptr;
    MTL::ComputeCommandEncoder *computeCommandEncoder = nullptr;

    ~MetalResources() {
        if ( computeCommandEncoder ) computeCommandEncoder->release();
        if ( computeCommandBuffer ) computeCommandBuffer->release();
        if ( computePrimaryPipelineState ) computePrimaryPipelineState->release();
        if ( kernelFunction ) kernelFunction->release();
        if ( library ) library->release();
        if ( libraryPath ) libraryPath->release();
        if ( commandQueuePrimary ) commandQueuePrimary->release();
        if ( device ) device->release();
    }

    bool initialise( const std::string& kernelFunctionName ) {
        device = MTL::CreateSystemDefaultDevice();
        if (!device)
        {
            std::cerr << "Failed to create a device\n";
            return false;
        }

        commandQueuePrimary = device->newCommandQueue();
        if (!commandQueuePrimary)
        {
            std::cerr << "Failed to create command queue\n";
            return false;
        }

        libraryPath = NS::String::string(METAL_SHADER_METALLIB_PATH, NS::UTF8StringEncoding);
        library = device->newLibrary(libraryPath, &error);
        if (!library)
        {
            std::cerr << "Failed to load the library from path: " << METAL_SHADER_METALLIB_PATH << "\n";
            return false;
        }

        kernelFunction = library->newFunction(NS::String::string(kernelFunctionName.c_str(), NS::UTF8StringEncoding));
        if (!kernelFunction)
        {
            std::cerr << "Failed to retrieve kernel function: " << kernelFunctionName << " Error: "
                      << error->localizedDescription()->utf8String() << "\n";
            return false;
        }

        computePrimaryPipelineState = device->newComputePipelineState(kernelFunction, &error);
        if (!computePrimaryPipelineState)
        {
            std::cerr << "Failed to set compute pipeline state. Error: "
                      << error->localizedDescription()->utf8String() << "\n";
            return false;
        }

        return true;
    }

    bool initialise( std::string const &kernelFunctionName, SimulationParameters const &keySimParams, SimulationStates const &keySimStates, SimulationFlags const &keySimFlags ) {

        device = MTL::CreateSystemDefaultDevice();
        if (!device)
        {
            std::cerr << "Failed to create a device\n";
            return false;
        }

        commandQueuePrimary = device->newCommandQueue();
        if (!commandQueuePrimary)
        {
            std::cerr << "Failed to create command queue\n";
            return false;
        }

        libraryPath = NS::String::string(METAL_SHADER_METALLIB_PATH, NS::UTF8StringEncoding);
        library = device->newLibrary(libraryPath, &error);
        if (!library)
        {
            std::cerr << "Failed to load the library from path: " << METAL_SHADER_METALLIB_PATH << "\n";
            return false;
        }

        // Define function constants here
        auto *constantValues = MTL::FunctionConstantValues::alloc()->init();
        constantValues->setConstantValue(&keySimParams, MTL::DataTypeStruct, static_cast<NS::Integer>(0));
        constantValues->setConstantValue(&keySimStates, MTL::DataTypeStruct, static_cast<NS::Integer>(1));
        constantValues->setConstantValue(&keySimFlags, MTL::DataTypeStruct, static_cast<NS::Integer>(2));

        auto *kernelFunctionString = NS::String::string(kernelFunctionName.c_str(), NS::UTF8StringEncoding);
        kernelFunction = library->newFunction(kernelFunctionString, constantValues, &error);
        if (!kernelFunction)
        {
            std::cerr << "Failed to retrieve kernel function: " << kernelFunctionName << " Error: "
                      << error->localizedDescription()->utf8String() << "\n";
            return false;
        }

        computePrimaryPipelineState = device->newComputePipelineState(kernelFunction, &error);
        if (!computePrimaryPipelineState)
        {
            std::cerr << "Failed to set compute pipeline state. Error: "
                      << error->localizedDescription()->utf8String() << "\n";
            return false;
        }

        return true;
    }

    bool manualRelease() {
        if ( computeCommandEncoder ) computeCommandEncoder->release();
        if ( computeCommandBuffer ) computeCommandBuffer->release();
        if ( computePrimaryPipelineState ) computePrimaryPipelineState->release();
        if ( kernelFunction ) kernelFunction->release();
        if ( library ) library->release();
        if ( libraryPath ) libraryPath->release();
        if ( commandQueuePrimary ) commandQueuePrimary->release();
        if ( device ) device->release();
        return true;
    }
};

#endif //SPINCHAINS_LOCALMETALRESOURCES_H
