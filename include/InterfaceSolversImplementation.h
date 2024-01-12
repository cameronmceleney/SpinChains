//
// Created by Cameron Aidan McEleney on 19/11/2023.
//

#ifndef SPINCHAINS_INTERFACESOLVERSIMPLEMENTATION_H
#define SPINCHAINS_INTERFACESOLVERSIMPLEMENTATION_H

class SimulationParameters;

class iSolversImplementation {
public:
    virtual void runMethod() = 0;
    virtual void                _testShockwaveConditions(double iteration) = 0;
            // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    virtual void                SolveRK2Classic() = 0;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    virtual void                SolveRK2() = 0;

    virtual void                RK2Parallel() = 0;

    virtual void                RK4Parallel() = 0;

    ~iSolversImplementation() = default;

};


#endif //SPINCHAINS_INTERFACESOLVERSIMPLEMENTATION_H
