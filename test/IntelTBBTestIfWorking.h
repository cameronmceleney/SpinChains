//
// Created by Cameron McEleney on 21/11/2023.
//

#ifndef SPINCHAINS_INTELTBBTESTIFWORKING_H
#define SPINCHAINS_INTELTBBTESTIFWORKING_H



// C++ Standard Libraries
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

// C++ Third Party Libraries
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>


class IntelTBBTestIfWorking {
private:
    double _sineSumParallel(const std::vector<double> &data);

    double _sineSumSeries(const std::vector<double> &data);

public:
    void initialiseTesting();

};
#endif //SPINCHAINS_INTELTBBTESTIFWORKING_H
