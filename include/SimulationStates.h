//
// Created by Cameron Aidan McEleney on 20/11/2023.
//

#ifndef SPINCHAINS_SIMULATIONSTATES_H
#define SPINCHAINS_SIMULATIONSTATES_H

// C++ Standard Libraries
#include <vector>

class SimulationStates {
public:
    std::vector<double> gilbertVector = {0};
    std::vector<std::vector<double>> gilbertVectorMulti = {};
    std::vector<double> largestMNormMulti = {1e-50, 1e-50};
    std::vector<double> mx0 = {0};
    std::vector<double> my0 = {0};
    std::vector<double> mz0 = {0};

    // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins
    std::vector<double> exchangeVec;



    // Description missing
    std::vector<int> layerSpinPairs;

    // Description missing
    std::vector<int>    layerSpinsInChain;

    // Description missing
    std::vector<int>    layerTotalSpins;

    // Vectors containing magnetic components (m), along each axis, at the initial conditions for all spins. Leave as zero!

    std::vector<std::vector<std::vector<double>>> m0Nest;
    std::vector<std::vector<std::vector<double>>> m1Nest;
    std::vector<std::vector<std::vector<double>>> m2Nest;
};

#endif //SPINCHAINS_SIMULATIONSTATES_H
