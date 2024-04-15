//
// Created by Cameron Aidan McEleney on 03/04/2024.
//

#ifndef SPINCHAINS_COMMONTOOLS_H
#define SPINCHAINS_COMMONTOOLS_H

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

class CommonTools {
public:
    struct CoreParamsForDevice {
        size_t NUM_ELEMENTS;
        size_t NUM_ELEMENTS_WITH_DATA;
        size_t NUM_ELEMENTS_PADDING;
        ushort NUM_ELEMENTS_PER_THREAD;
        ushort NUM_VEC4;
        size_t bytes;

        // Method to calculate the size of the struct
        void calculateSize() {
            // Since the struct contains static-sized elements only,
            // we can directly use sizeof to get its total size.
            // This includes all its members, such as 'N' and 'bytes'.
            bytes = sizeof(*this);
        }

        // Constructor to automatically calculate bytes upon creation
        explicit CoreParamsForDevice(size_t N, size_t numDp) : NUM_ELEMENTS_WITH_DATA(N), NUM_ELEMENTS_PER_THREAD(numDp) {
            NUM_VEC4 = std::ceil(NUM_ELEMENTS_PER_THREAD / 4);  // ensures we have at least one float 4 per thread
            NUM_ELEMENTS_PADDING = 0;
            NUM_ELEMENTS = NUM_ELEMENTS_WITH_DATA;
            calculateSize();
        }

        explicit CoreParamsForDevice(size_t N_true, size_t N_pad,size_t numDp) : NUM_ELEMENTS_WITH_DATA(N_true),
        NUM_ELEMENTS_PADDING(N_pad), NUM_ELEMENTS_PER_THREAD(numDp) {
            NUM_VEC4 = std::ceil(NUM_ELEMENTS_PER_THREAD / 4);
            NUM_ELEMENTS = NUM_ELEMENTS_WITH_DATA + NUM_ELEMENTS_PADDING;
            calculateSize();
        }
    };

    struct CoreParamsForHost {
        size_t NUM_ELEMENTS = -1;
        size_t NUM_ELEMENTS_WITH_DATA = -1;
        size_t NUM_ELEMENTS_PADDING = -1;
        ushort NUM_THREADS = -1;
        size_t NUM_THREADGROUPS = -1;
        ushort DATAPOINTS_PER_CONTAINER = -1;
        ushort CONTAINERS_PER_THREAD = -1;
        ushort DATAPOINTS_PER_THREAD = -1;
        ushort SIZE_DTYPE = -1;
        size_t BYTES_FOR_ELEMENTS = -1;

        void calculateDataPerThread(size_t containersPerThread, size_t datapointsPerContainer)
        {
            CONTAINERS_PER_THREAD = containersPerThread;
            DATAPOINTS_PER_CONTAINER = datapointsPerContainer;
            DATAPOINTS_PER_THREAD = datapointsPerContainer * containersPerThread;
        }

        void calculateSize() {

            if (NUM_ELEMENTS_PADDING >= 0)
            {
                NUM_ELEMENTS = NUM_ELEMENTS_WITH_DATA + NUM_ELEMENTS_PADDING;
            }
            else
            {
                std::cout << "Padding less than 0 during host calculateSize." << std::endl;
                std::exit(1);
            }

            BYTES_FOR_ELEMENTS = NUM_ELEMENTS * SIZE_DTYPE;
        }

        // Constructor to automatically calculate bytes upon creation
        explicit CoreParamsForHost(size_t N) : NUM_ELEMENTS_WITH_DATA(N) {
            SIZE_DTYPE = sizeof(float);
            NUM_ELEMENTS_PADDING = 0;
            calculateDataPerThread(1, 1);
            calculateSize();
        }
        CoreParamsForHost() = default;
    };

public:
    static std::vector<float> getRandomVector(size_t size);

    static void padVectorForCacheLine(std::vector<float>& vec, CoreParamsForHost &coreParams);

    static void verifyResultsAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                              const std::vector<float> &resultVec, const bool &printAnswer);

    static void verifyResultsSineAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                              const std::vector<float> &resultVec, const bool &printAnswer);

    static void verifyResultsSineAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                          const std::vector<float> &resultVec, CoreParamsForDevice &coreParams, const bool &printAnswer);

    static void verifyResultsSineAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                      const float *resultVec, CoreParamsForDevice &coreParams, const bool &printAnswer);

};


#endif //SPINCHAINS_COMMONTOOLS_H
