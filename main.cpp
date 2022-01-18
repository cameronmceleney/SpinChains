#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"
#include "GlobalVariables.h"

int main() {

    GlobalVariables GV;

    SpinChainEigenSolverClass SolverClass;
    Numerical_Methods_Class NMMethods;

    bool query = false;

    int in_numSpins; // number of spins in the chain
    std::cout << "Enter the number of spins in the chain: ";
    std::cin >> in_numSpins; // Takes user input for the number of spins
    GV.SetNumSpins(in_numSpins);

    double in_exchangeMin;
    std::cout << "Enter the minimum exchange value: ";
    std::cin >> in_exchangeMin;
    GV.SetExchangeMinVal(in_exchangeMin);
    std::cout << "Main minval: " << GV.GetExchangeMinVal() << std::endl;

    double in_exchangeMax;
    std::cout << "Enter the maximum exchange value: ";
    std::cin >> in_exchangeMax;
    GV.SetExchangeMaxVal(in_exchangeMax);
    std::cout << "Main maxval: " << GV.GetExchangeMaxVal() << std::endl;

    std::string in_fileNameBase; //Better name might be fileID
    std::cout << "Enter the unique identifier that all filenames will share: ";
    std::cin >> in_fileNameBase;
    GV.SetFileNameBase(in_fileNameBase);

    if (query) {
        std::string in_filePath;
        std::cout << "Enter the absolute path to the directory to save data: ";
        std::cin >> in_filePath;
        GV.SetFilePath(in_filePath);
    }

    if (query) {
        SolverClass.SolveInputs();
    }
    NMMethods.NMSetup();
    NMMethods.RK2();

    return 0;

}