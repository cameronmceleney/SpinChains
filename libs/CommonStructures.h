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

    // Axis definitions for improved readability in indexing
    struct Axis {
        static constexpr std::size_t x = 0;
        static constexpr std::size_t y = 1;
        static constexpr std::size_t z = 2;
    };


    // VectorND: Generalized N-dimensional vector with direct member access
    template <typename T, std::size_t N>
    struct VectorND {
        static_assert(N > 0 && N <= 3, "VectorND only supports dimensions 1, 2, or 3");

        std::array<T, N> data;

        // Constructor initializes array with the provided arguments
        template <typename... Args, typename = std::enable_if_t<sizeof...(Args) == N>>
        constexpr VectorND(Args... args) : data{{args...}} {}

        // Direct member access
        constexpr T& x() {
            static_assert(N >= 1, "Access to member 'x' is out of bounds for this vector size");
            return data[0];
        }

        constexpr const T& x() const {
            static_assert(N >= 1, "Access to member 'x' is out of bounds for this vector size");
            return data[0];
        }

        constexpr T& y() {
            static_assert(N >= 2, "Access to member 'y' is out of bounds for this vector size");
            return data[1];
        }

        constexpr const T& y() const {
            static_assert(N >= 2, "Access to member 'y' is out of bounds for this vector size");
            return data[1];
        }

        constexpr T& z() {
            static_assert(N >= 3, "Access to member 'z' is out of bounds for this vector size");
            return data[2];
        }

        constexpr const T& z() const {
            static_assert(N >= 3, "Access to member 'z' is out of bounds for this vector size");
            return data[2];
        }

        // Access operators
        constexpr T& operator[](std::size_t index) {
            return data[index];
        }

        constexpr const T& operator[](std::size_t index) const {
            return data[index];
        }

        constexpr bool operator==(const VectorND& other) const {
            return data == other.data;
        }

        constexpr bool operator<(const VectorND& other) const {
            return data < other.data;
        }
    };

    // PointND: Represents a point in N-dimensional space, leveraging VectorND for position and value
    template <typename PositionType, std::size_t PDim, typename ValueType, std::size_t VDim>
    struct PointND {
        VectorND<PositionType, PDim> pos;  // position
        VectorND<ValueType, VDim> val;     // value

        // Constructor using initializer lists
        constexpr PointND(VectorND<PositionType, PDim> position, VectorND<ValueType, VDim> value)
            : pos(position), val(value) {}

        // Constructor using individual arguments for both position and value
        template <typename... PosArgs, typename... ValArgs>
        constexpr PointND(PosArgs... posArgs, ValArgs... valArgs)
            : pos{posArgs...}, val{valArgs...} {
            static_assert(sizeof...(PosArgs) == PDim, "Incorrect number of arguments for position");
            static_assert(sizeof...(ValArgs) == VDim, "Incorrect number of arguments for value");
        }

        // Overload operator[] to access value via index
        constexpr ValueType& operator[](std::size_t axis) {
            return val[axis];
        }

        constexpr const ValueType& operator[](std::size_t axis) const {
            return val[axis];
        }
    };

    // Timer utility for measuring elapsed time
    struct Timer {
        std::string timerName;
        std::chrono::time_point<std::chrono::system_clock> startSolver;
        std::chrono::time_point<std::chrono::system_clock> endSolver;
        size_t solverElapsedTime = 0;
        bool useMilliseconds = false;

        void start(const bool &useMs = false) {
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

        void setName(const std::string &name) {
            timerName = name;
        }

        void cleanUp() {
            solverElapsedTime = 0;
            useMilliseconds = false;
        }

        void print() {
            auto printTimePoint = [](auto &timePoint) {
                std::time_t printTime = std::chrono::system_clock::to_time_t(timePoint);
                return std::ctime(&printTime);
            };

            std::cout << "----------------------------------------------------------------\n";
            if (!timerName.empty())
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

    // Units of measurement
    enum class Unit {
        Tesla,
        Gauss,
        JoulesPerMeter
    };

    // Parallelization methods
    enum class Parallelisations {
        Sequential,
        Multithreaded
    };

    // Type aliases for common point configurations
    using Point1D = PointND<int, 1, double, 1>;
    using Point2D = PointND<int, 2, double, 2>;
    using Point3D = PointND<int, 3, double, 3>;

    using Vector1D = VectorND<double, 1>;
    using Vector2D = VectorND<double, 2>;
    using Vector3D = VectorND<double, 3>;

}

// Provide easier access to Axis without polluting the global namespace
using Axis = CommonStructures::Axis;

#endif //SPINCHAINS_COMMONSTRUCTURES_H
