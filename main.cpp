#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"
#include "CommonLibs.h"


int main() {

    // GitHub Token: ***REMOVED*** (works as of 14 Oct 22)
    SpinChainEigenSolverClass SolverClass{};
    Numerical_Methods_Class RK2_method_use{};

    // Select between eigenvalue derivation and numerical modelling
    bool findEigenvalues = true;
    GV.SetEmailWhenCompleted(false);

    // Set global file-related parameters
    GV.SetCurrentTime();

    // Set global simulation parameters
    GV.SetAnisotropyField(0.0);
    GV.SetStaticBiasField(1e-3);
    GV.SetNumSpins(300);
    GV.SetExchangeMinVal(21);
    GV.SetExchangeMaxVal(21);
    GV.SetGyromagneticConstant(28.2);
    GV.SetIsFerromagnetic(true);

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
    GV.SetFilePath("MacOS", findEigenvalues);

    if (findEigenvalues) {
        // std::cout << "Finding eigenvalues and eigenvectors" << std::endl;
        SolverClass.CalculateEigFreqs();
    }

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