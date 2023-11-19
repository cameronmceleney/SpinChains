//
// Created by Cameron McEleney on 31/10/2023.
//

#ifndef SPINCHAINS_NMDATAHANDLING_H
#define SPINCHAINS_NMDATAHANDLING_H

// C++ User Libraries (Parent Class)
#include "../include/NMSuperClassTest.h"
#include "NMMethods.h"

class NMMethods;

class NMDataHandling:  public NMSuperClassTest {
    friend class NMMethods;
public:
    NMDataHandling(std::shared_ptr<SystemDataContainer> data) : NMSuperClassTest(data){};
    void                performInitialisation() override { std::cout<<"Data Handling called correctly" << std::endl; };

public:
    // ####################################            Define Private Variables            ###################################

    // ####################################            Define Private Variables            ###################################
    // Description missing
    std::vector<double> flattenNestedVector(const std::vector<std::vector<double>>& nestedVector);

    // ####################################    Functions to control data output    ####################################
    // Description missing
    void                CreateColumnHeaders(std::ofstream &outputFileName);

    // Description missing
    void                CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata=false);

    // Description missing
    void                InformUserOfCodeType(const std::string& nameNumericalMethod);

    // Description missing
    void                PrintVector(std::vector<double> &vectorToPrint, bool shouldExitAfterPrint);

    // Description missing
    void                PrintNestedNestedVector(std::vector<std::vector<std::vector<double>>> nestedNestedVector);

    // Description missing
    void                SaveDataToFile(std::ofstream &outputFileName, std::vector<double> &arrayToWrite, int &iteration);

    // Description missing
    void                SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration);

    // Description missing
    void                CreateMetadata(bool print_end_time=false);
        /*
     * ################################################################################################################
     * #########    Overloaded functions to control data output that are specific to multilayered systems    ##########
     * ################################################################################################################
     */

    // Description missing
    void CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata, int layer);

    // Description missing
    void SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration, int layer);

    // Description missing
    void CreateColumnHeaders(std::ofstream &outputFileName, int& layer);

//protected:
    // ####################################            Define Protected Variables            ###################################

    // ####################################            Define Protected Functions            ###################################

//public:
    // ####################################            Define Public Variables            ###################################

    // ####################################            Define Public Functions            ###################################

};


#endif //SPINCHAINS_NMDATAHANDLING_H
