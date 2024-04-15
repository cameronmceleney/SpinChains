#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>
#include <string>

class Timer {
public:
    Timer();
    void start(bool useMs = false);
    void stop();
    void setName(const std::string& name);
    void print() const;
    void cleanUp();

private:
    std::string timerName;
    std::chrono::time_point<std::chrono::high_resolution_clock> startSolver, endSolver;
    long long solverElapsedTime;
    bool useMilliseconds;
};

Timer::Timer() : solverElapsedTime(0), useMilliseconds(false) {}

void Timer::start(bool useMs) {
    useMilliseconds = useMs;
    startSolver = std::chrono::high_resolution_clock::now();
}

void Timer::stop() {
    endSolver = std::chrono::high_resolution_clock::now();
    solverElapsedTime = useMilliseconds ?
        std::chrono::duration_cast<std::chrono::milliseconds>(endSolver - startSolver).count() :
        std::chrono::duration_cast<std::chrono::seconds>(endSolver - startSolver).count();
}

void Timer::setName(const std::string& name) {
    timerName = name;
}

void Timer::print() const {
    std::cout << "----------------------------------------------------------------\n";
    if (!timerName.empty())
        std::cout << "Timing Information - " << timerName << "\n";
    else
        std::cout << "\nTiming Information. \n";

    std::cout << "\tElapsed: " << solverElapsedTime;
    if (useMilliseconds)
        std::cout << " [milliseconds]\n";
    else
        std::cout << " [seconds]\n";

    std::cout << "----------------------------------------------------------------\n";
}

void Timer::cleanUp() {
    timerName.clear();
    solverElapsedTime = 0;
    useMilliseconds = false;
}

#endif // TIMER_H