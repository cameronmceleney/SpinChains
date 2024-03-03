//
// Created by Cameron Aidan McEleney on 20/11/2023.
//

#ifndef SPINCHAINS_SIMULATIONSTATES_H
#define SPINCHAINS_SIMULATIONSTATES_H

// C++ Standard Libraries
#include <list>
#include <map>
#include <vector>

class SimulationStates {
public:
    /**
     * Holds a linearly space array of values to describe all the exchange damping factors
     */
    std::vector<double> gilbertVector = {0};

    /**
     * Description missing
     */
    std::vector<std::vector<double>> gilbertVectorMulti = {};

    /**
     * Description missing
     */
    std::vector<double> largestMNormMulti = {1e-50, 1e-50};

    /**
     * Magnetic moments along x-axis (m_x) for all spins at the initial stage of each timestep
     */
    std::vector<double> mx0 = {0};

    /**
     * Magnetic moments along y-axis (m_y) for all spins at the initial stage of each timestep
     */
    std::vector<double> my0 = {0};

    /**
     * Magnetic moments along z-axis (m_z) for all spins at the initial stage of each timestep
     */
    std::vector<double> mz0 = {0};

    /**
     * Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins
     */
    std::vector<double> exchangeVec;

    /**
     * Sites to be printed if SimulationParameters->shouldPrintDiscreteSites is true
     */
    std::list <int> fixedOutputSites;

    /**
     * Description missing
     */
    std::vector<int> layerSpinPairs;

    /**
     * Description missing
     */
    std::vector<int>    layerSpinsInChain;

    /**
     * Description missing
     */
    std::vector<int>    layerTotalSpins;

    std::vector<std::vector<std::vector<double>>> m0Nest;
    std::vector<std::vector<std::vector<double>>> m1Nest;
    std::vector<std::vector<std::vector<double>>> m2Nest;

    /**
     * Description missing
     */
     std::vector<int> discreteDrivenSites;

     std::map<int, std::pair<double, int>> dRGradientMap;
     std::map<int, std::pair<double, int>> dmiGradientMap;
     std::map<int, std::pair<double, int>> dampingGradientMap;
};

#endif //SPINCHAINS_SIMULATIONSTATES_H
