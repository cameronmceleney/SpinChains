//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_SOLVERSIMPLEMENTATION_H
#define SPINCHAINS_SOLVERSIMPLEMENTATION_H

// C++ Corresponding Interface
#include "InterfaceSolversImplementation.h"

// C++ Standard Libraries
#include <atomic>
#include <string>
#include <vector>

// C++ Third Party Library
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

// C++ User Libraries (General)
#include "GlobalVariables.h"
#include "../libs/CommonStructures.h"
#include "../libs/progressbar.hpp"

// C++ User Libraries (Parent Class)
#include "SolversSuperClass.h"
// C++ User Libraries (Sibling Classes)
#include "SolversDataHandling.h"

// C++ User Libraries (Class' Components)
#include "DemagnetisationFields.h"
#include "DzyaloshinskiiMoriyaInteraction.h"
#include "EffectiveFields.h"
#include "MagnetisationDynamics.h"
#include "DipolarFields.h"
#include "SpinTransferTorque.h"
#include "ExchangeField.h"
#include "BiasFields.h"

class SolversImplementation :
        public SolversSuperClass, public iSolversImplementation {
private:
    DemagnetisationFields demagField;
    DzyaloshinskiiMoriyaInteraction dmInteraction;
    EffectiveFields effectiveField;
    MagnetisationDynamics llg;
    DipolarFields dipolarField;
    SpinTransferTorque stt;
    ExchangeField exchangeField;
    BiasFields biasField;

private:
    int progressBarSubdivisions = 100;
    void _resizeClassContainers();
    void _resizeClassContainersTest();
    void _resizeClassContainersRK4();

    static void _testOutputValues( double& mxTerm, double& myTerm, double& mzTerm, int site, int iteration, const std::string &rkStage );
    /**
     * Description missing
     */
    void _testShockwaveConditions( double iteration ) override;

    void    _clearDataStorageVariables();

    /**
     * Transfer results to a regular vector outside the parallel region
     */
    static void _transferDataThenReleaseAtomicVector( std::vector<std::atomic<double>> &atomicVector,
                                                      std::vector<double> &regularVector,
                                                      bool shouldRelease = true);

    /**
     * Description missing
     */
    void _clearThenReleaseVector( auto &vec );

    /**
     * Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method. Original
     * implementation of the RK2 method; used for testing purposes
     */
    void SolveRK2Classic() override;

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void SolveRK2() override;

    /**
     * Evaluate the given system using the RK2 method which has been configured to enabled
     * parallelisation through Intel's oneAPI TBB
     */
    void RK2Parallel() override;

    /**
     * Evaluate the given system using the RK2 method which has been configured to enabled
     * parallelisation through Intel's oneAPI TBB
     */
    void RK4Parallel() override;

    /**
     * Description missing
     */
    void RK2StageMultithreaded( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                std::vector<double> &myOut, std::vector<double> &mzOut,
                                double &currentTime, double &stepsize, int &iteration,
                                std::string rkStage );
    void RK2StageMultithreadedTest( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                    const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                    std::vector<double> &myOut, std::vector<double> &mzOut,
                                    const double &currentTime, const double &stepsize,
                                    const int &iteration, std::string rkStage );

    void RK2StageMultithreadedCompact( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                    const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                    std::vector<double> &myOut, std::vector<double> &mzOut,
                                    const double &currentTime, const double &stepsize,
                                    const int &iteration, std::string rkStage );

    void RK4StageMultithreadedCompact( const std::vector<double> &mxIn, const std::vector<double> &myIn,
                                    const std::vector<double> &mzIn, std::vector<double> &mxOut,
                                    std::vector<double> &myOut, std::vector<double> &mzOut,
                                    const double &currentTime, const double &stepsize,
                                    const int &iteration, std::string rkStage );

public:
    SolversImplementation( std::shared_ptr<SimulationManager> sharedSimManager,
                           std::shared_ptr<SimulationParameters> paramsData,
                           std::shared_ptr<SimulationStates> sharedSimStates,
                           std::shared_ptr<SimulationFlags> sharedSimFlags );

    ~SolversImplementation() = default;
public:
    // TODO. Temp as public for testing. In future, will be private members (need to sort inheritance first)
    // TODO. Update all methods (currently only Parallel) to use class members as containers
    std::vector<double> effectiveFieldX;
    std::vector<double> effectiveFieldY;
    std::vector<double> effectiveFieldZ;

    std::vector<double> demagXp;
    std::vector<double> demagYp;
    std::vector<double> demagZp;
    std::vector<double> dipoleXp;
    std::vector<double> dipoleYp;
    std::vector<double> dipoleZp;
    std::vector<double> dmiXp;
    std::vector<double> dmiYp;
    std::vector<double> dmiZp;
    std::vector<double> sttXp;
    std::vector<double> sttYp;
    std::vector<double> sttZp;

    // RK Stage Containers for reuse
    std::vector<double> mx1p;
    std::vector<double> my1p;
    std::vector<double> mz1p;
    std::vector<double> mx2p;
    std::vector<double> my2p;
    std::vector<double> mz2p;
    std::vector<double> mx3p;
    std::vector<double> my3p;
    std::vector<double> mz3p;
    std::vector<double> mx4p;
    std::vector<double> my4p;
    std::vector<double> mz4p;

    std::vector<double> mxk;
    std::vector<double> myk;
    std::vector<double> mzk;
public:
    struct DipoleTerms {
        double x;
        double y;
        double z;

        DipoleTerms() : x(0.0), y(0.0), z(0.0){}
        DipoleTerms(double updateX, double updateY, double updateZ) : x(updateX), y(updateY), z(updateZ){}
    };

    struct DMITerms {
        double x;
        double y;
        double z;

        DMITerms() : x(0.0), y(0.0), z(0.0){}
        DMITerms(double updateX, double updateY, double updateZ) : x(updateX), y(updateY), z(updateZ){}
    };

    struct DemagTerms {
        double x;
        double y;
        double z;

        DemagTerms() : x(0.0), y(0.0), z(0.0){}
        DemagTerms(double updateX, double updateY, double updateZ) : x(updateX), y(updateY), z(updateZ){}
    };

    struct STTTerms {
        double x;
        double y;
        double z;

        STTTerms() : x(0.0), y(0.0), z(0.0){}
        STTTerms(double updateX, double updateY, double updateZ) : x(updateX), y(updateY), z(updateZ){}
    };

    struct HkTerms {
        double x;
        double y;
        double z;

        HkTerms() : x(0.0), y(0.0), z(0.0){}
        HkTerms(double updateX, double updateY, double updateZ) : x(updateX), y(updateY), z(updateZ){}
    };

    struct MkTerms {
        double x;
        double y;
        double z;

        MkTerms() : x(0.0), y(0.0), z(0.0){}
        MkTerms(double updateX, double updateY, double updateZ) : x(updateX), y(updateY), z(updateZ){}
    };
public:
    void performInitialisation() override { runMethod(); };
    void reinitialise() override { runMethod(); };

    void runMethod() override;
};


#endif //SPINCHAINS_SOLVERSIMPLEMENTATION_H
