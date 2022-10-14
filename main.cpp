#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"
#include "CommonLibs.h"

int main() {

    // GitHub Token: ***REMOVED*** (works as of 14 Oct 22)
    SpinChainEigenSolverClass SolverClass{};
    Numerical_Methods_Class RK2_method_use{};

    // Set core parameters
    GV.SetCurrentTime();
    GV.SetFilePath();

    // Set simulation parameters
    GV.SetAnisotropyField(0.787);
    GV.SetStaticBiasField(0.1);
    GV.SetNumSpins(static_cast<int>(500));
    GV.SetExchangeMinVal(53);
    GV.SetExchangeMaxVal(53);
    GV.SetGyromagneticConstant(28.8E9);

    // I keep forgetting to check the exchanges, hence this warning
    if (GV.GetExchangeMinVal() == GV.GetExchangeMaxVal()) {
        std::cout << "Uniform Exchange\n";
    } else {
        std::cout << "Non-Uniform Exchange\n";
    }

    /* int in_numSpins; // number of spins in the chain
     * std::cout << "Enter the number of spins in the chain: ";
     * std::cin >> in_numSpins; // Takes user input for the number of spins
     * GV.SetNumSpins(in_numSpins);

     * double in_exchangeMin;
     * std::cout << "Enter the minimum exchange value: ";
     * std::cin >> in_exchangeMin;
     * GV.SetExchangeMinVal(in_exchangeMin)

     * double in_exchangeMax;
     * std::cout << "Enter the maximum exchange value: ";
     * std::cin >> in_exchangeMax;
     * GV.SetExchangeMaxVal(in_exchangeMax);
     */

    std::string in_fileNameBase; //Better name might be fileID
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> in_fileNameBase;
    GV.SetFileNameBase("T"+in_fileNameBase);

    /* std::string in_filePath;
     * std::cout << "Enter the absolute path to the directory to save data: ";
     * std::cin >> in_filePath;
     * GV.SetFilePath(in_filePath);
     */

    // SolverClass.CalculateEigFreqs();
    RK2_method_use.NMSetup();
    RK2_method_use.RK2MidpointAFM();

    return 0;
}