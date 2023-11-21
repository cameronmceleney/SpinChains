
// Corresponding header
#include "../include/GlobalVariables.h"

double GlobalVariablesClass::GetAnisotropyField() {
    return _anisotropyField;
}
void GlobalVariablesClass::SetAnisotropyField(double anisotropyField) {
    _anisotropyField = anisotropyField;
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

bool GlobalVariablesClass::GetEmailWhenCompleted() {
    if (_shouldSendEmail.has_value()) {
        if (_shouldSendEmail.value()) {
            std::cout << "This functionality is not yet implemented" << std::endl;
            std::exit(1);
        } else
            return _shouldSendEmail.value();
    } else {
        std::cout << "Boolean value for _shouldSendEmail not set" << std::endl;
        std::exit(1);
    }
}
void GlobalVariablesClass::SetEmailWhenCompleted(bool shouldSendEmail) {
    _shouldSendEmail = shouldSendEmail;
}

double GlobalVariablesClass::GetExchangeMaxVal() {
    return _exchangeMaxVal;
}
void GlobalVariablesClass::SetExchangeMaxVal(double exchangeMaxVal) {
    _exchangeMaxVal = exchangeMaxVal;
}

double GlobalVariablesClass::GetExchangeMinVal() {
    return _exchangeMinVal;
}
void GlobalVariablesClass::SetExchangeMinVal(double exchangeMinVal) {
    _exchangeMinVal = exchangeMinVal;
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

std::string GlobalVariablesClass::GetFileNameBase() {
    return _fileNameBase;
}
void GlobalVariablesClass::SetFileNameBase(std::string fileNameBase) {
    _fileNameBase = fileNameBase;
}

std::string GlobalVariablesClass::GetFilePath() {
    return _filePath;
}
void GlobalVariablesClass::SetFilePath(std::string osName) {

    for (size_t i = 0; i < osName.size(); ++i) {
        osName[i] = toupper(static_cast<unsigned char>(osName[i]));
    }

    if (GetShouldFindEigenvalues()) {
        std::filesystem::path filepath = _filePath;
        bool filepathExists = std::filesystem::is_directory(filepath.parent_path());

        if (!filepathExists) {
            // Guard clause
            std::cout << "The filepath " << filepath << " does not exist. Please check the filepath and try again." << std::endl;
        }

        std::string filenameExtension = _fileNameBase + "_Eigens";
        std::filesystem::create_directory(_filePath + filenameExtension);
        _filePath += filenameExtension + "/";
    } else {
        if (osName == "MACOS") {
            // Default Windows filepath for my laptop
            _filePath = "/Users/cameronaidanmceleney/CLionProjects/Data/" + FindDateToday() + "/Simulation_Data/";
        } else if (osName == "WINDOWS") {
            // Default Windows filepath for my desktop
            _filePath = "D:/Data/" + FindDateToday() + "/Simulation_Data/";
        } else {
            // Guard clause
            std::cout << "The operating system name " << osName << " is not recognised. Please check the operating system name and try again." << std::endl;
            std::exit(1);
        }
    }
}

double GlobalVariablesClass::GetGyromagneticConstant() {
    return _gyromagneticConstant;
}
void GlobalVariablesClass::SetGyromagneticConstant(double gyromagneticConstant) {
    _gyromagneticConstant = gyromagneticConstant * 1e9 * 2 * M_PI;
}

bool GlobalVariablesClass::GetIsFerromagnetic() {
    return _isFerromagnetic;
}
void GlobalVariablesClass::SetIsFerromagnetic(bool isFerromagnetic) {
    _isFerromagnetic = isFerromagnetic;
}

bool GlobalVariablesClass::GetIsExchangeUniform() {
    return _isExchangeUniform;
}
void GlobalVariablesClass::SetIsExchangeUniform() {
    if (_exchangeMaxVal == -3.141592 || _exchangeMinVal == -3.141592) {
        std::cout << "One or more exchange values are not set. Please check the values and try again." << std::endl;
        std::exit(1);
    } else if (_exchangeMaxVal == _exchangeMinVal)
        _isExchangeUniform = true;
    else
        _isExchangeUniform = false;
}

int GlobalVariablesClass::GetNumSpins() {
    return _numSpins;
}
void GlobalVariablesClass::SetNumSpins(int numSpins) {
    _numSpins = numSpins;
}

double GlobalVariablesClass::GetStaticBiasField() {
    return _staticBiasField;
}
void GlobalVariablesClass::SetStaticBiasField(double staticBiasField) {
    _staticBiasField = staticBiasField;
}

std::string GlobalVariablesClass::GetNumericalMethod() {
    return _chosenNumericalMethod;
}
void GlobalVariablesClass::SetNumericalMethod(std::string chosenNumericalMethod) {
    _chosenNumericalMethod = chosenNumericalMethod;
}

bool GlobalVariablesClass::GetShouldFindEigenvalues() {
    if (_shouldFindEigenvalues.has_value())
        return _shouldFindEigenvalues.value();
    else {
        std::cout << "Boolean value for _shouldFindEigenvalues not set" << std::endl;
        std::exit(1);
    }
}
void GlobalVariablesClass::SetShouldFindEigenvalues(bool shouldFindEigenvalues) {
    _shouldFindEigenvalues = shouldFindEigenvalues;
}