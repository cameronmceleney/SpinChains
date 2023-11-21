//
// Created by Cameron Aidan McEleney on 19/11/2023.
//

#ifndef SPINCHAINS_INMMETHODS_H
#define SPINCHAINS_INMMETHODS_H

#include <memory>

class SimulationParameters;

class INMMethods {
public:
    virtual void runMethod() = 0;
    virtual void                _testShockwaveConditions(double iteration) = 0;
            // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    virtual void                SolveRK2Classic() = 0;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    virtual void                SolveRK2() = 0;

    ~INMMethods() = default;

};


#endif //SPINCHAINS_INMMETHODS_H
