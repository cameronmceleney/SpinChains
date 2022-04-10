#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.h"
#include "CommonLibs.h"

int main() {

    // GitHub Token: ***REMOVED*** (works as of 22 Mar 22)
    SpinChainEigenSolverClass SolverClass{};
    Numerical_Methods_Class RK2_method_use{};

    GV.SetFilePath();
    GV.SetBiasField(0.1);
    GV.SetNumSpins(static_cast<int>(3));
    // GV.SetNumSpins(static_cast<int>(6000 * (1.0 + 0.05 * 2)));
    GV.SetExchangeMinVal(43.5);
    GV.SetExchangeMaxVal(43.5);


    /* int in_numSpins; // number of spins in the chain
     * std::cout << "Enter the number of spins in the chain: ";
     * std::cin >> in_numSpins; // Takes user input for the number of spins
     * GV.SetNumSpins(in_numSpins);
     */

    /* double in_exchangeMin;
     * std::cout << "Enter the minimum exchange value: ";
     * std::cin >> in_exchangeMin;
     * GV.SetExchangeMinVal(in_exchangeMin)
     */
    /* double in_exchangeMax;
     * std::cout << "Enter the maximum exchange value: ";
     * std::cin >> in_exchangeMax;
     * GV.SetExchangeMaxVal(in_exchangeMax);
     */

    std::string in_fileNameBase; //Better name might be fileID
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> in_fileNameBase;
    GV.SetFileNameBase("LLGTest"+in_fileNameBase);

    /* std::string in_filePath;
     * std::cout << "Enter the absolute path to the directory to save data: ";
     * std::cin >> in_filePath;
     * GV.SetFilePath(in_filePath);
     */

    /* SolverClass.CalculateEigFreqs();
     * SolverClass.CalculateEigFreqs();
     */

    RK2_method_use.NMSetup();
    RK2_method_use.RK2LLG();

    return 0;
}