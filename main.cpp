#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"
#include "CommonLibs.h"

int main() {

    // GitHub Token: ***REMOVED*** (works as of 14 Oct 22)
    SpinChainEigenSolverClass SolverClass{};
    Numerical_Methods_Class RK2_method_use{};

    // Set global file-related parameters
    GV.SetCurrentTime();
    GV.SetFilePath();

    // Set global simulation parameters
    GV.SetAnisotropyField(0);
    GV.SetStaticBiasField(0.1);
    GV.SetNumSpins(4000);
    GV.SetExchangeMinVal(43.5);
    GV.SetExchangeMaxVal(132.0);
    GV.SetGyromagneticConstant(29.2e9);
    GV.SetIsFerromagnetic(true);

    // Select between eigenvalue derivation and numerical modelling
    bool findEigenvalues = false;

    // I keep forgetting to check the exchanges, hence this warning
    if (GV.GetExchangeMinVal() == GV.GetExchangeMaxVal()) {
        std::cout << "Uniform Exchange\n";
    } else {
        std::cout << "Non-Uniform Exchange\n";
    }

    std::string in_fileNameBase; //Better name might be fileID
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> in_fileNameBase;
    GV.SetFileNameBase("T"+in_fileNameBase);

    RK2_method_use.NMSetup();
    RK2_method_use.RK2MidpointFM();
    exit(0);

    if (findEigenvalues)
        SolverClass.CalculateEigFreqs();

#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnreachableCode"
    else if (!findEigenvalues) {
        RK2_method_use.NMSetup();
        if (GV.GetIsFerromagnetic())
            RK2_method_use.RK2MidpointFM();
        else
            RK2_method_use.RK2MidpointAFM();
    }
#pragma clang diagnostic pop
    return 0;
}