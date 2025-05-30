#ifndef SPINCHAINS_GLOBALVARIABLES_H
#define SPINCHAINS_GLOBALVARIABLES_H

// C++ Standard Libraries
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>  // Use to create sentinels for booleans
#include <string>

class GlobalVariablesClass{

private:
//  Dtype               Member Name                             // Comment
    double              _anisotropyField;                       // Anisotropy field [T]
    double              _dmiConstant;

    std::string         _currentTime;
    std::string         _dateToday;                             // Contains today's date.
    double              _exchangeMaxVal = -3.141592;            // Maximum exchange field constant value (J_min).
    double              _exchangeMinVal = -3.141592;            // minimum exchange field constant value (J_min).
    bool                _isExchangeUniform;                     //

    std::string         _fileNameBase;                          // Forms the first part of any filename generated by one iteration of this software.
    std::string         _filePath;                              // Path to directory where all files will be saved.
    double              _gyromagneticConstant;                  // Gyromagnetic constant [2 Pi GHz/T].
    bool                _isFerromagnetic;                       // If the material is ferromagnetic [true] or antiferromagnetic [false]
    int                 _numSpins;                              // Number of sites (spins) in the chain.

    std::optional<bool>  _shouldFindEigenvalues;
    std::optional<bool>  _shouldSendEmail;
    double              _staticBiasField;                       // Applied bias field [T]. Often written as H_0 or H_static in literature.

    std::string         _chosenNumericalMethod;
    // Private Functions
    std::string         FindDateToday();                        // Finds today's date in the standard format for this code.

public:
//  Dtype               Member Name                                     // Comment
    double              GetAnisotropyField();                           // Get the anisotropy field [T].
    void                SetAnisotropyField(double anisotropyField);     // Set the anisotropy field [T]. Often called H_A or H_a.

    int                 GetNumSpins();                                  // Get the total number of spins in chain.
    void                SetNumSpins(int numSpins);                      // Set the number of spins in chain. Must always be updated every time a new region is added.

    std::string         GetFilePath();                                  // Get file path to today's directory.
    /**
     * Set file path to directory for writing data.
     *
     * Date is automatically found, but the remaining directory tree must be manually inserted.
     *
     * @param os_name String variable that should be either "MacOS" or "Windows".
     * @param isEigenValues Boolean variable to indicate if the file path is for eigenvalues or numerical modelling.
     */
    void                SetFilePath(std::string osName);

    std::string         GetFileNameBase();                              // Get the custom suffix (base) name.
    void                SetFileNameBase(std::string fileNameBase);      // Set a custom suffix to all filenames for this simulation.

    bool                GetEmailWhenCompleted();
    void                SetEmailWhenCompleted(bool shouldSendEmail);
    double              GetExchangeMinVal();                            // Get the minimum exchange interaction strength [T].
    void                SetExchangeMinVal(double exchangeMinimum);      // Set the minimum exchange interaction strength [T].

    double              GetExchangeMaxVal();                            // Get the maximum exchange interaction strength [T].
    void                SetExchangeMaxVal(double exchangeMaximum);      // Set the maximum exchange interaction strength H_ex [T].

    double              GetStaticBiasField();                           // Get the applied external field [T].
    void                SetStaticBiasField(double staticBiasField);     // Set the applied external field [T]. Often called H_0 or H_{ext}.

    std::string         GetCurrentTime();                               // Get the time at the moment this method was set.
    void                SetCurrentTime();                               // Finds the current time (ISO8601 format) at the moment the method was invoked.

    double              GetGyromagneticConstant();                      // Get the gyromagnetic ratio
    void                SetGyromagneticConstant(double gyromagneticConstant);   // Set the gyromagnetic ratio, called 'gamma'. Enter in [GHz] as conversion is handled by method!

    bool                GetIsFerromagnetic();                           // Get the material type
    void                SetIsFerromagnetic(bool isFerromagnetic);       // Set if the material is: ferromagnetic [true] or antiferromagnetic [false]

    bool                GetIsExchangeUniform();                         // Get if the exchange interaction is uniform or not.
    void                SetIsExchangeUniform();                         // Set if the exchange interaction is uniform or not.

    bool                GetShouldFindEigenvalues();
    void                SetShouldFindEigenvalues(bool shouldFindEigenvalues);

    std::string         GetNumericalMethod();

    /**
     * Set the numerical method to be used for the simulation.
     * @param chosenNumericalMethod Can be either "RK2", "RK2c", or "RK2p"
     */
    void                SetNumericalMethod(std::string chosenNumericalMethod);

    double GetDMIConstant();
    void SetDMIConstant( double dmiConstant );

    std::string getUserName();

};


#endif //SPINCHAINS_GLOBALVARIABLES_H
