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
void GlobalVariablesClass::SetFilePath(std::string filePath) {
    _filePath = std::move(filePath);
}

std::string GlobalVariablesClass::GetFileNameBase() {
    return _fileNameBase;
}
void GlobalVariablesClass::SetFileNameBase(std::string fileNameBase) {
    _fileNameBase = std::move(fileNameBase);
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
