#include "GlobalVariables.h"

int GlobalVariablesClass::GetNumSpins() {
    return _numSpins;
}
void GlobalVariablesClass::SetNumSpins(int numSpins) {
    _numSpins = numSpins;
}

std::string GlobalVariablesClass::GetFilePath() {
    return _filePath;
}
void GlobalVariablesClass::SetFilePath() {
    _filePath = "/Users/cameronmceleney/CLionProjects/Data/"+ FindDateToday() +"/RK2 Shockwaves Tests Data/"; // This filepath is for Mac!
    //_filePath = "D:/Data/" + FindDateToday() +"/RK2 Shockwaves Tests Data/" // This filepath is for Windows
}

std::string GlobalVariablesClass::GetFileNameBase() {
    return _fileNameBase;
}
void GlobalVariablesClass::SetFileNameBase(std::string fileNameBase) {
    _fileNameBase = fileNameBase;
}

double GlobalVariablesClass::GetExchangeMinVal() {
    return _exchangeMinVal;
}
void GlobalVariablesClass::SetExchangeMinVal(double exchangeMinVal) {
    _exchangeMinVal = exchangeMinVal;
}

double GlobalVariablesClass::GetExchangeMaxVal() {
    return _exchangeMaxVal;
}
void GlobalVariablesClass::SetExchangeMaxVal(double exchangeMaxVal) {
    _exchangeMaxVal = exchangeMaxVal;
}

std::string GlobalVariablesClass::FindDateToday() {

    long t = std::time(nullptr);
    tm tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%d %b %y");
    std::string date = oss.str();

    return date;
}
