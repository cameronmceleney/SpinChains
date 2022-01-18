#include "GlobalVariables.h"

int GlobalVariables::GetNumSpins() {
    return _numSpins;
}
void GlobalVariables::SetNumSpins(int numSpins) {
    _numSpins = numSpins;
}

std::string GlobalVariables::GetFilePath() {
    return _filePath;
}
void GlobalVariables::SetFilePath(std::string filePath) {
    _filePath = filePath;
}

std::string GlobalVariables::GetFileNameBase() {
    return _fileNameBase;
}
void GlobalVariables::SetFileNameBase(std::string fileNameBase) {
    _fileNameBase = fileNameBase;
}

double GlobalVariables::GetExchangeMinVal() {
    return _exchangeMinVal;
}
void GlobalVariables::SetExchangeMinVal(double exchangeMinVal) {
    _exchangeMinVal = exchangeMinVal;
}

double GlobalVariables::GetExchangeMaxVal() {
    return _exchangeMaxVal;
}
void GlobalVariables::SetExchangeMaxVal(double exchangeMaxVal) {
    _exchangeMaxVal = exchangeMaxVal;
}
