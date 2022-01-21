#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"
#include "CommonLibs.h"

int main() {

    SpinChainEigenSolverClass SolverClass{};
    Numerical_Methods_Class RK2_method_use{};

    bool query = false;

    //int in_numSpins; // number of spins in the chain
    //std::cout << "Enter the number of spins in the chain: ";
    //std::cin >> in_numSpins; // Takes user input for the number of spins
    //GV.SetNumSpins(in_numSpins);
    GV.SetNumSpins(4000);
    std::cout << "Values are hardcoded as: Jmin=70T, Jmax=70T, nspins=4000, DRStart=1, stepsize=1e-15, itermax=3.5e5."<<std::endl;

    //double in_exchangeMin;
    //std::cout << "Enter the minimum exchange value: ";
    //std::cin >> in_exchangeMin;
    //GV.SetExchangeMinVal(in_exchangeMin);
    GV.SetExchangeMinVal(70);


    //double in_exchangeMax;
    //std::cout << "Enter the maximum exchange value: ";
    //std::cin >> in_exchangeMax;
    //GV.SetExchangeMaxVal(in_exchangeMax);
    GV.SetExchangeMaxVal(70);


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
        SolverClass.CalculateEigFreqs();
    }
    RK2_method_use.NMSetup();
    RK2_method_use.RK2Shockwaves();

    return 0;

}