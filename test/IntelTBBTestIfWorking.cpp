//
// Created by Cameron McEleney on 21/11/2023.
//

// Corresponding Header
#include "IntelTBBTestIfWorking.h"

double IntelTBBTestIfWorking::_sineSumParallel(const std::vector<double>& data) {
    return tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0, data.size()),
                    0.0,
                    [&data](const tbb::blocked_range<size_t>& r, double running_total) {
                        for (size_t i = r.begin(); i != r.end(); i++) {
                            running_total += std::sin(data[i]);
                        }
                        return running_total; // Combines the totals from different threads
                },
                std::plus<double>() // Reduce the results from all threads
    );
}

double IntelTBBTestIfWorking::_sineSumSeries(const std::vector<double> &data) {
    double running_total = 0.0;
    for (size_t i = 0; i < data.size(); ++i)
        running_total += std::sin(data[i]);

    return running_total;
}

void IntelTBBTestIfWorking::initialiseTesting() {
    // Create a large vector of values to then compute the sign of
    std::vector<double> dataToProcess(4000);
    for (size_t i = 0; i < dataToProcess.size(); ++i)
        dataToProcess[i] = static_cast<double>(i);
    std::vector<double> dataSeries = dataToProcess;

    // Start timing for parallel calculation
    auto startParallel = std::chrono::high_resolution_clock::now();

    // Call the parallel sine sum method
    double parallelResult = _sineSumParallel(dataToProcess);

    // Stop timing for parallel calculation
    auto stopParallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedParallel = stopParallel - startParallel;
    std::cout << "Parallel Sine Sum (" << parallelResult << ") took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(elapsedParallel).count() << "ms" << std::endl;

    // Start timing for serial calculation
    auto startSeries = std::chrono::high_resolution_clock::now();

    // Call the series sine sum method
    double seriesResult = _sineSumSeries(dataSeries);

    // Stop timing for serial calculation
    auto stopSeries = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSeries = stopSeries - startSeries;
    std::cout << "Series Sine Sum (" << seriesResult << ") took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(elapsedSeries).count() << "ms" << std::endl;

    // Print the speedup
    std::cout << "Speedup: " << (elapsedSeries.count() / elapsedParallel.count())*100 << "%" << std::endl;
}