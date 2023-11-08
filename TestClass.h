//
// Created by Cameron McEleney on 08/11/2023.
//

#ifndef SPINCHAINS_TESTCLASS_H
#define SPINCHAINS_TESTCLASS_H

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>

class TestClass {
private:
    double  _parallelSineSum(const std::vector<double>& data);
    double  _parallelSineSumSeries(const std::vector<double>& data);

public:
    void    testFunction();
};


#endif //SPINCHAINS_TESTCLASS_H
