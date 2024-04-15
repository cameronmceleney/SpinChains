//
// Created by Cameron Aidan McEleney on 03/04/2024.
//

#include "CommonTools.h"

std::vector<float> CommonTools::getRandomVector(size_t size) {
    std::vector<float> vec(size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1.0);

    for (auto& el : vec) {
        el = static_cast<float>(dis(gen));
    }

    return vec;
}

void CommonTools::padVectorForCacheLine(std::vector<float>& vec, CoreParamsForHost &coreParams) {
    // ONLY WORKS FOR FLOAT FOR NOW DUE TO FUNCTION PROTOTYPE, AND PUSH_BACK(0.0F)

    size_t bytesCacheLine = 16;  // value in bytes
    // e.g. for a float bytesPerElement=4 so elementsPerCacheLine = 4
    size_t elementsPerCacheLine = bytesCacheLine / coreParams.SIZE_DTYPE;

    // Check number of elements in vector is a multiple of the cache line width
    size_t remainder = coreParams.NUM_ELEMENTS_WITH_DATA % elementsPerCacheLine;
    if (remainder != 0) {
        coreParams.NUM_ELEMENTS_PADDING = elementsPerCacheLine - remainder;
        for (size_t i = 0; i < coreParams.NUM_ELEMENTS_PADDING; ++i) {
            vec.push_back(0.0f); // Pad with zeros for now for a FLOAT
        }
    }
}

void CommonTools::verifyResultsAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                                  const std::vector<float> &resultVec, const bool &printAnswer) {

    for (size_t i = 0; i < inVecA.size(); ++i) {
        if (resultVec[i] != inVecA[i] + inVecB[i]) {
            std::cerr << "Verification failed at index " << i << ": " << resultVec[i] << " != " << inVecA[i] << " + "
            << inVecB[i] << std::endl;
            throw std::runtime_error("Result verification failed.");
        }

        if (printAnswer)
        {
            std::cout << resultVec[i] << " = " << inVecA[i] << " + " << inVecB[i] << std::endl;
        }
    }
}

void CommonTools::verifyResultsSineAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                                  const std::vector<float> &resultVec, const bool &printAnswer) {

    constexpr float epsilon = 0.000001;
    for (size_t i = 0; i < inVecA.size(); ++i) {
        float resultCpu = sin(inVecA[i] * inVecB[i]) + inVecA[i];
        float diffCpuToGpu = resultCpu - resultVec[i];
        if (fabs(diffCpuToGpu) >= epsilon) {
            // Beyond acceptable range
            std::cerr << "Verification failed at index " << i << " with diff " << diffCpuToGpu << " | " << resultVec[i] << " != sin(" << inVecA[i] << " * "
            << inVecB[i] << ") + " << inVecA[i] << std::endl;
            throw std::runtime_error("Result verification failed.");
        }

        if (printAnswer)
        {
            std::cout << resultVec[i] << " = " << inVecA[i] << " + " << inVecB[i] << std::endl;
        }
    }
}

void CommonTools::verifyResultsSineAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                                  const std::vector<float> &resultVec, CoreParamsForDevice &coreParams, const bool &printAnswer) {

    constexpr float epsilon = 0.000001;
    for (size_t i = 0; i < coreParams.NUM_ELEMENTS; ++i) {
        float resultCpu = sin(inVecA[i] * inVecB[i]) + inVecA[i];
        float diffCpuToGpu = resultCpu - resultVec[i];
        if (fabs(diffCpuToGpu) >= epsilon) {
            // Beyond acceptable range
            std::cerr << "Verification failed at index " << i << " with diff " << diffCpuToGpu << " | " << resultVec[i] << " != sin(" << inVecA[i] << " * "
            << inVecB[i] << ") + " << inVecA[i] << std::endl;
            throw std::runtime_error("Result verification failed.");
        }

        if (printAnswer)
        {
            std::cout << resultVec[i] << " = " << inVecA[i] << " + " << inVecB[i] << std::endl;
        }
    }
}

void CommonTools::verifyResultsSineAddition(const std::vector<float> &inVecA, const std::vector<float> &inVecB,
                                  const float *resultVec, CoreParamsForDevice &coreParams, const bool &printAnswer) {

    constexpr float epsilon = 0.000001;
    for (size_t i = 0; i < coreParams.NUM_ELEMENTS; ++i) {
        float resultCpu = sin(inVecA[i] * inVecB[i]) + inVecA[i];
        float diffCpuToGpu = resultCpu - resultVec[i];
        if (fabs(diffCpuToGpu) >= epsilon) {
            // Beyond acceptable range
            std::cerr << "Verification failed at index " << i << " with diff " << diffCpuToGpu << " | " << resultVec[i] << " != sin(" << inVecA[i] << " * "
            << inVecB[i] << ") + " << inVecA[i] << std::endl;
            throw std::runtime_error("Result verification failed.");
        }

        if (printAnswer)
        {
            std::cout << resultVec[i] << " = " << inVecA[i] << " + " << inVecB[i] << std::endl;
        }
    }
}