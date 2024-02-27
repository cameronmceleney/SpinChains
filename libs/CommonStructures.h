//
// Created by Cameron McEleney on 27/02/2024.
//

#ifndef SPINCHAINS_COMMONSTRUCTURES_H
#define SPINCHAINS_COMMONSTRUCTURES_H

// C++ Standard Libraries
#include <array>
#include <chrono>
#include <iostream>

namespace CommonStructures {
    // Name stands for "User Common Structures"

    struct Point3D {
        std::array<int, 3> coordinates;

        explicit constexpr Point3D(int x, int y, int z) : coordinates{{x, y, z}} {}

        bool operator==(const Point3D& other) const {
            return coordinates == other.coordinates;
        }

        bool operator<(const Point3D& other) const {
            return coordinates < other.coordinates;
        }

        // Convenience accessors
        int x() const { return coordinates[0]; }
        int y() const { return coordinates[1]; }
        int z() const { return coordinates[2]; }
    };

    struct Timer {
        std::string timerName;
        std::chrono::time_point<std::chrono::system_clock> startSolver;
        std::chrono::time_point<std::chrono::system_clock> endSolver;
        long long solverElapsedTime = 0;
        bool useMilliseconds = false;

        void start(bool useMs = false) {
            if (useMs)
                useMilliseconds = true;
            startSolver = std::chrono::system_clock::now();
        }

        void stop() {
            endSolver = std::chrono::system_clock::now();
            if (useMilliseconds)
                solverElapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endSolver - startSolver).count();
            else
                solverElapsedTime = std::chrono::duration_cast<std::chrono::seconds>(endSolver - startSolver).count();
        }

        void setName( const std::string &name ) {
            timerName = name;
        }

        void cleanUp() {
            solverElapsedTime = 0;
            useMilliseconds = false;
        }

        void print() {
            auto printTimePoint = []( auto &timePoint) {
                                        std::time_t printTime = std::chrono::system_clock::to_time_t(timePoint);
                                        return std::ctime(&printTime);
            };

            std::cout << "----------------------------------------------------------------\n";
            if ( !timerName.empty())
                std::cout << "Timing Information - " << timerName << "\n";
            else
                std::cout << "\nTiming Information. \n";

            std::cout << "\tStart: " << printTimePoint(startSolver);
            std::cout << "\tEnd: " << printTimePoint(endSolver);
            if (useMilliseconds)
                std::cout << "\tElapsed: " << solverElapsedTime << " [milliseconds]" << std::endl;
            else
                std::cout << "\tElapsed: " << solverElapsedTime << " [seconds]" << std::endl;

            std::cout << "----------------------------------------------------------------\n";

            cleanUp();
        }
    };
}

#endif //SPINCHAINS_COMMONSTRUCTURES_H
