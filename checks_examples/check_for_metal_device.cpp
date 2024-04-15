//
// check_for_metal_device.cpp
// Created by Cameron Aidan McEleney on 08/03/2024.
//

#include "check_for_metal_device.h"
#include "CustomTimer.hpp"
#include "LocalMetalResources.h"

#include <iostream>
#include <chrono>
#include <iostream>
#include <random>


class CoffeeExample::CoffeeExampleImpl {
private:

    std::vector<float> vec1;
    std::vector<float> vec2;

    bool FULL_DEBUG;

public:
    // Here you would define all the methods and member variables that were previously in CommonTools.
    // For simplicity, I'll only demonstrate a few.

    void vectorAddition( const int &numElements, const std::string &kernelFunctionName,
                                        const bool &fullDebug) {
        CoreParamsForHost hostCoreParams(numElements);
        hostCoreParams.calculateDataPerThread(1, 1);

        FULL_DEBUG = fullDebug;
        Timer gpuTimer, kernelTimer;
        gpuTimer.setName("Overall Timer (Vector Addition)(Shared Resources)");
        kernelTimer.setName("Kernel Timer (Vector Addition)(Shared Resources)");

        // Original code
        //auto vec1 = getRandomVector(hostCoreParams.NUM_ELEMENTS);
        //auto vec2 = getRandomVector(hostCoreParams.NUM_ELEMENTS);

        // Initialise random numbers in each array
        if (vec1.empty() || vec1.size() != hostCoreParams.NUM_ELEMENTS_WITH_DATA)
            vec1 = getRandomVector(hostCoreParams.NUM_ELEMENTS_WITH_DATA);
        if  (vec2.empty() || vec2.size() != hostCoreParams.NUM_ELEMENTS_WITH_DATA)
            vec2 = getRandomVector(hostCoreParams.NUM_ELEMENTS_WITH_DATA);

        std::vector<float> vec3(hostCoreParams.NUM_ELEMENTS_WITH_DATA);

        // Exclude costly vector population as this isn't part of the timings tests.
        gpuTimer.start(true);

        CoreParamsForDevice deviceCoreParams(hostCoreParams.NUM_ELEMENTS_WITH_DATA, hostCoreParams.DATAPOINTS_PER_THREAD);

        if (FULL_DEBUG)
        {
            std::cout << "Elements in vectors. (vec1): " << vec1.size() << " | (vec2): " << vec2.size() << std::endl;
        }

        // Initialise Metal device, command queue, and pipeline
        MetalResources localMTLRes;
        localMTLRes.initialise(kernelFunctionName);

        // Allocate shared memory for the host and device
        auto bufferVec1 = localMTLRes.device->newBuffer(hostCoreParams.BYTES_FOR_ELEMENTS, MTL::ResourceStorageModeShared);
        auto bufferVec2 = localMTLRes.device->newBuffer(hostCoreParams.BYTES_FOR_ELEMENTS, MTL::ResourceStorageModeShared);
        auto bufferVec3 = localMTLRes.device->newBuffer(hostCoreParams.BYTES_FOR_ELEMENTS, MTL::ResourceStorageModeShared);
        auto bufferCoreParameters = localMTLRes.device->newBuffer(deviceCoreParams.bytes, MTL::ResourceStorageModeShared);

        // Error handling example for buffer creation
        if ( !bufferVec1 || !bufferVec2 || !bufferVec3 || !bufferCoreParameters )
        {
            std::cerr << "Failed to create one or more buffers!" << std::endl;
            std::exit(1);
        }

        // Encode commands
        localMTLRes.computeCommandBuffer = localMTLRes.commandQueuePrimary->commandBuffer();
        localMTLRes.computeCommandEncoder = localMTLRes.computeCommandBuffer->computeCommandEncoder();
        localMTLRes.computeCommandEncoder->setComputePipelineState(localMTLRes.computePrimaryPipelineState);

        localMTLRes.computeCommandEncoder->setBuffer(bufferVec1, 0, 0);
        localMTLRes.computeCommandEncoder->setBuffer(bufferVec2, 0, 1);
        localMTLRes.computeCommandEncoder->setBuffer(bufferVec3, 0, 2);
        localMTLRes.computeCommandEncoder->setBuffer(bufferCoreParameters, 0, 3);

        // Copy data from host to device (CPU -> GPU)
        memcpy(bufferVec1->contents(), vec1.data(), hostCoreParams.BYTES_FOR_ELEMENTS);
        memcpy(bufferVec2->contents(), vec2.data(), hostCoreParams.BYTES_FOR_ELEMENTS);
        memcpy(bufferCoreParameters->contents(), &deviceCoreParams, deviceCoreParams.bytes);

        // Threads per block (1<<10 == 1024)
        hostCoreParams.NUM_THREADS = 1 << 10;
        if ( FULL_DEBUG && (hostCoreParams.NUM_THREADS > localMTLRes.device->maxThreadsPerThreadgroup().height ) )
        {
            std::cout << "Attempted to set more threads per thread-group than device permits\n";
            std::exit(1);
        }
        MTL::Size NUM_THREADS_PER_THREADGROUP = {hostCoreParams.NUM_THREADS, 1, 1};

        // Blocks per Grid (padded as required)
        hostCoreParams.NUM_THREADGROUPS = (hostCoreParams.NUM_ELEMENTS + hostCoreParams.NUM_THREADS - 1) / hostCoreParams.NUM_THREADS;
        MTL::Size NUM_THREADGROUPS_PER_GRID = {hostCoreParams.NUM_THREADGROUPS, 1, 1};

        kernelTimer.start(true);
        // Inform the kernel of the layout of the threads (similar to <<<NUM_BLOCKS, NUM_THREADS>>> in CUDA
        localMTLRes.computeCommandEncoder->dispatchThreadgroups(NUM_THREADGROUPS_PER_GRID, NUM_THREADS_PER_THREADGROUP);
        // Send instruction that our pipeline state encoding is completed
        localMTLRes.computeCommandEncoder->endEncoding();

        // Asynchronously launch the kernel on the device (GPU)
        //commandBuffer->GPUStartTime()
        std::cout << "\n\n";
        localMTLRes.computeCommandBuffer->commit();

        // Ensure synchronisation is explicitly managed (don't want to rely on memcpy for now)
        localMTLRes.computeCommandBuffer->waitUntilCompleted();
        kernelTimer.stop();
        kernelTimer.print();

        // Copy sum vector from device to host (GPU -> CPU)
        memcpy(vec3.data(), bufferVec3->contents(), hostCoreParams.BYTES_FOR_ELEMENTS);

        // No need to time checking results; not part of test
        gpuTimer.stop();
        gpuTimer.print();

        // Check result for errors
        verifyResultsSineAddition(vec1, vec2, vec3, false);

        // Free memory on the device (GPU)
        bufferVec1->release();
        bufferVec2->release();
        bufferVec3->release();

        // Release all other resources
        localMTLRes.manualRelease();
    }

    // Add additional methods and member variables as needed.
};

// Constructor
CoffeeExample::CoffeeExample() : pCoffeeExampleImpl(std::make_unique<CoffeeExampleImpl>()) {}

// Destructor needs to be defined where Impl is a complete type, which is here in the .cpp file.
CoffeeExample::~CoffeeExample() = default;

// Public method implementations that delegate to the pImpl
void CoffeeExample::vectorAddition( const int &numElements, const std::string &kernelFunctionName,
                                    const bool &fullDebug ) {
    pCoffeeExampleImpl->vectorAddition(numElements, kernelFunctionName, fullDebug);
}

// Implement other methods similarly, delegating to the Impl class.