//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSSUPERCLASS_H
#define SPINCHAINS_SOLVERSSUPERCLASS_H

// C++ Standard Libraries
//#include <algorithm>  // Will need in future. Left here so I don't forget

// C++ Third Party Libraries

// C++ User Libraries (General)
#include "CommonLibs.h"

// C++ User Libraries (Containers)
#include "SimulationFlags.h"
#include "SimulationParameters.h"
#include "SimulationStates.h"

class SolversSuperClass {
public:
    struct CustomTimer {
        std::string timerName;
        std::chrono::time_point<std::chrono::system_clock> startSolver;
        std::chrono::time_point<std::chrono::system_clock> endSolver;
        long long solverElapsedTime = 0;

        void start() {
            startSolver = std::chrono::system_clock::now();
        }

        void stop() {
            endSolver = std::chrono::system_clock::now();
            solverElapsedTime = std::chrono::duration_cast<std::chrono::seconds>(endSolver - startSolver).count();
        }

        void setName( const std::string &name ) {
            timerName = name;
        }

        void print() const {
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
            std::cout << "\tElapsed: " << solverElapsedTime << " [seconds]" << std::endl;
            std::cout << "----------------------------------------------------------------\n";
        }
    };

protected:
    // Getter for child classes to access sharedVariables
    std::shared_ptr<SimulationParameters> simParams;
    std::shared_ptr<SimulationStates> simStates;
    std::shared_ptr<SimulationFlags> simFlags;
    CustomTimer methodTimer;


public:
    //SolversSuperClass(std::shared_ptr<SimulationParameters> simParams);
    SolversSuperClass( std::shared_ptr<SimulationParameters> sharedSimParams,
                       std::shared_ptr<SimulationStates> sharedSimStates,
                       std::shared_ptr<SimulationFlags> sharedSimFlags );

    virtual ~SolversSuperClass(); // Destructor (though not strictly necessary with smart pointers)

public:
    virtual void performInitialisation() = 0;

    static std::shared_ptr<SolversSuperClass>
    createSimulationInstance( std::shared_ptr<SimulationParameters> sharedSimParams,
                              std::shared_ptr<SimulationStates> sharedSimStates,
                              std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SolversSuperClass>
    createConfigurationInstance( std::shared_ptr<SimulationParameters> sharedSimParams,
                                 std::shared_ptr<SimulationStates> sharedSimStates,
                                 std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SolversSuperClass>
    createDataHandlingInstance( std::shared_ptr<SimulationParameters> sharedSimParams,
                                std::shared_ptr<SimulationStates> sharedSimStates,
                                std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SolversSuperClass>
    createMethodsInstance( std::shared_ptr<SimulationParameters> sharedSimParams,
                           std::shared_ptr<SimulationStates> sharedSimStates,
                           std::shared_ptr<SimulationFlags> sharedSimFlags );

    static std::shared_ptr<SimulationParameters> getSimulationParamsContainer();

    static std::shared_ptr<SimulationStates> getSimulationStatesContainer();

    static std::shared_ptr<SimulationFlags> getSimulationFlagsContainer();

    static std::pair<double, double> testInstance();
};


#endif //SPINCHAINS_SOLVERSSUPERCLASS_H
