//
// Created by Cameron McEleney on 08/11/2023.
//

#include "TestClass.h"

double TestClass::_parallelSineSum(const std::vector<double>& data) {
    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, data.size()),
        0.0,
        [&data](const tbb::blocked_range<size_t>& r, double running_total) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                running_total += std::sin(data[i]);
            }
            return running_total; // Combine totals from different threads.
        },
        std::plus<double>() // Reduce the results from all threads.
    );
}

double TestClass::_parallelSineSumSeries(const std::vector<double>& data) {
    double running_total = 0.0;

    // Serial iteration over the vector
    for (size_t i = 0; i < data.size(); ++i) {
        running_total += std::sin(data[i]);
    }

    return running_total; // Return the total sum
}



void TestClass::testFunction() {
    // Create a large vector of values to calculate the sine of.
    std::vector<double> dataParallel(4000);
    for (size_t i = 0; i < dataParallel.size(); ++i) {
        dataParallel[i] = static_cast<double>(i);
    }
    std::vector<double> dataSeries = dataParallel;

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Call the parallel sine sum function.
    double sumParallel = _parallelSineSum(dataParallel);

    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedParallel = end - start;
    std::cout << "The sum of sines is: " << sumParallel << std::endl;
    std::cout << "Parallel calculation took: " << elapsedParallel.count() << " seconds." << std::endl;

    // Start timing for serial calculation
    start = std::chrono::high_resolution_clock::now();

    // Call the serial sine sum function.
    double sumSeries = _parallelSineSumSeries(dataSeries);

    // Stop timing
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSeries = end - start;
    std::cout << "The sum of sines is: " << sumSeries << std::endl;
    std::cout << "Serial calculation took: " << elapsedSeries.count() << " seconds." << std::endl;
    std::cout << "Speedup: " << elapsedSeries.count() / elapsedParallel.count() << std::endl;
}