#include "GlobalVariables.h"

double GlobalVariablesClass::GetBiasField() {
    return _biasField;
}
void GlobalVariablesClass::SetBiasField(double biasField) {
    _biasField = biasField;
}

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
    _filePath = "/Users/cameronmceleney/CLionProjects/Data/"+ FindDateToday() +"/Simulation_Data/"; // This filepath is for Mac!
    // _filePath = "D:/Data/" + FindDateToday() +"/Simulation_Data/"; // This filepath is for Windows
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

std::string GlobalVariablesClass::GetCurrentTime() {
    return _currentTime;
}
void GlobalVariablesClass::SetCurrentTime() {
    time_t     now = time(0);
    struct tm  timeStruct;
    char       buf[80];
    timeStruct = *localtime(&now);

    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format

    strftime(buf, sizeof(buf), "T%H%M", &timeStruct);

    _currentTime = buf;
}

std::string GlobalVariablesClass::FindDateToday() {

    time_t     now = time(0);
    struct tm  timeStruct;
    char       buf[80];
    timeStruct = *localtime(&now);

    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format

    strftime(buf, sizeof(buf), "%Y-%m-%d", &timeStruct);

    return buf;
}

