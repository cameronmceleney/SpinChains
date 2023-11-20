//
// Created by Cameron McEleney on 31/10/2023.
//

#include "NMDataHandling.h"

void NMDataHandling::CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata) {
    /**
     * Write all non-data information to the output file.
     */
    if (is_metadata) {
        outputFileName << "Key Data\n\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using magDynamics: [" << systemData->useLLG << "]\t\t\t\tUsing Shockwave: [" << systemData->hasShockwave << "]\t\tDrive from LHS: [" << systemData->lhsDrive <<
                       "]\nNumerical Method Used: [" << methodUsed << "]\t\tHas Static Drive: [" << systemData->hasStaticDrive << "]\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0): " << GV.GetStaticBiasField() << " T\t\t\t" << "Dynamic Bias Field (H_D1): " << systemData->dynamicBiasField << " T\n" <<
                          "Dynamic Bias Field Scale Factor: " << systemData->shockwaveInitialStrength << "\t\t" << "Second Dynamic Bias Field (H_D2): " << systemData->shockwaveMax << " T\n" <<
                          "Driving Frequency (f): " << systemData->drivingFreq << "Hz\t\t""Driving Region Start Site: " << systemData->drivingRegionLhs - systemData->numSpinsDamped << "\n" <<
                          "Driving Region End Site: " << systemData->drivingRegionRhs - systemData->numSpinsDamped << " \t\t\t" << "Driving Region Width: " << systemData->drivingRegionWidth << " \n" <<
                          "Max. Sim. Time: " << systemData->maxSimTime << " s\t\t\t\t" << "Min. Exchange Val (J): " << systemData->exchangeEnergyMin  << " T\n" <<
                          "Max. Exchange Val (J): " << systemData->exchangeEnergyMax << " T\t\t\t" << "Max. Iterations: " << systemData->iterationEnd << "\n" <<
                          "No. DataPoints: " << systemData->numberOfDataPoints << " \t\t\t\t" << "No. Spins in Chain: " << systemData->numSpinsInChain << "\n" <<
                          "No. Damped Spins: " << systemData->numSpinsDamped << "per side\t\t\t" << "No. Total Spins: " << systemData->systemTotalSpins << " \n" <<
                          "systemData->stepsize (h): " << systemData->stepsize << "\t\t\t\t" << "Gilbert Damping Factor: " << systemData->gilbertDamping << "\n" <<
                          "Gyromagnetic Ratio (2Pi*Y): " << systemData->gyroMagConst << "\t\t""Shockwave Gradient Time: " << systemData->iterStartShock << "s\n" <<
                          "Shockwave Application Time: " << systemData->shockwaveGradientTime * systemData->stepsize << "s\n" <<
                          std::endl;

        return;
    }
    else {

        outputFileName << "Key Data\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using magDynamics," << systemData->useLLG << ",Using Shockwave," << systemData->hasShockwave << ",Drive from LHS," << systemData->lhsDrive <<
                       ",Numerical Method Used," << methodUsed << ",Has Static Drive," << systemData->hasStaticDrive << "\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                          "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site, Driving Region Width,"
                          "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                          "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, systemData->stepsize (h),Gilbert Damping Factor, Gyromagnetic Ratio (2Pi*Y),"
                          "Shockwave Gradient Time [s], Shockwave Application Time [s]"
                          "\n";

        outputFileName << GV.GetStaticBiasField() << ", " << systemData->dynamicBiasField << ", " << systemData->shockwaveInitialStrength << ", " << systemData->shockwaveMax << ", "
                       << systemData->drivingFreq << ", " << systemData->drivingRegionLhs - systemData->numSpinsDamped << ", " << systemData->drivingRegionRhs - systemData->numSpinsDamped << ", " << systemData->drivingRegionWidth << ", "
                       << systemData->maxSimTime << ", " << systemData->exchangeEnergyMin << ", " << systemData->exchangeEnergyMax << ", " << systemData->iterationEnd << ", " << systemData->numberOfDataPoints << ", "
                       << systemData->numSpinsInChain << ", " << systemData->numSpinsDamped << ", " << systemData->systemTotalSpins << ", " << systemData->stepsize << ", " << systemData->gilbertDamping << ", " << systemData->gyroMagConst << ", "
                       << systemData->iterStartShock << ", " << systemData->shockwaveGradientTime * systemData->stepsize
                       << "\n";

        outputFileName << "\n";
    }

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments );
    outputFileName << "Note(s):," << notesComments << "\n"; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n";

    outputFileName << "\n";

    CreateColumnHeaders(outputFileName);

    std::cout << "\n";
}
void NMDataHandling::CreateColumnHeaders(std::ofstream &outputFileName) {
    /**
     * Creates the column headers for each spin site simulated. This code can change often, so compartmentalising it in
     * a separate function is necessary to reduce bugs.
     */
    if (systemData->printAllData or systemData->printFixedLines) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for (int i = 1; i <= systemData->systemTotalSpins; i++) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if (systemData->printFixedSites) {

        outputFileName << "Time";
        for (int & fixed_out_val : systemData->fixedOutputSites)
            outputFileName << "," << fixed_out_val;
        outputFileName << std::endl;

        //outputFileName << "Time" << ", "
        //               << static_cast<int>(14000) << ","
        //               << static_cast<int>(16000) << ","
        //               << static_cast<int>(18000) << ","
        //               << static_cast<int>(20000) << std::endl;

    }
}
void NMDataHandling::InformUserOfCodeType(const std::string& nameNumericalMethod) {
    /**
     * Informs the user of the code type they are running, including: solver type; special modules.
     */
    if (systemData->useLLG)
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (magDynamics) code";
    else
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (Torque) code";

    if (systemData->hasShockwave)
        std::cout << " with shockwave module.\n";
    else
        std::cout << ".\n";

}
void NMDataHandling::PrintVector(std::vector<double> &vectorToPrint, bool shouldExitAfterPrint) {

    std::cout << "\n\n";

    int count = 0;
    for (double i: vectorToPrint) {
        if (++count % 10 == 0)
            std::cout << std::setw(8) << i << std::endl;
        else
            std::cout << std::setw(8) << i << ", ";
    }
    std::cout << "\n\n";

    if (shouldExitAfterPrint)
        exit(0);
}
void NMDataHandling::PrintNestedNestedVector(std::vector<std::vector<std::vector<double>>> nestedNestedVector){
        // Print the contents of nestedNestedVector1
    std::cout << "{";
    for (size_t i = 0; i < nestedNestedVector.size(); ++i) {
        std::cout << "{";
        for (size_t j = 0; j < nestedNestedVector[i].size(); ++j) {
            std::cout << "{";
            for (size_t k = 0; k < nestedNestedVector[i][j].size(); ++k) {
                std::cout << nestedNestedVector[i][j][k];
                if (k < nestedNestedVector[i][j].size() - 1) {
                    std::cout << ",";
                }
            }
            std::cout << "}";
            if (j < nestedNestedVector[i].size() - 1) {
                std::cout << ",";
            }
        }
        std::cout << "}";
        if (i < nestedNestedVector.size() - 1) {
            std::cout << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}
void NMDataHandling::SaveDataToFile(std::ofstream &outputFileName, std::vector<double> &arrayToWrite, int &iteration) {
    std::cout.precision(6);
    std::cout << std::scientific;

    if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
        if (systemData->printFixedLines) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * systemData->stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if (systemData->printFixedSites) {
            /*outputFileName << (iteration * systemData->stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * systemData->stepsize);
            for (int & fixed_out_val : systemData->fixedOutputSites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (systemData->printAllData) {
        for (int i = 0; i <= systemData->systemTotalSpins; i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * systemData->stepsize) << ","; // Print current time
            else if (i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (systemData->printFixedLines) {
        // iteration >= static_cast<int>(systemData->iterationEnd / 2.0) &&
        if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
            //if (iteration == systemData->iterationEnd) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * systemData->stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;
        }
    } else {
        if (systemData->printAllData) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * systemData->stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
                if (systemData->printFixedSites) {

                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[systemData->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[systemData->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[systemData->systemTotalSpins] << std::endl;

                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[systemData->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[systemData->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(systemData->systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(systemData->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[systemData->systemTotalSpins] << std::endl;
                }
            }
        }
    } */
}
void NMDataHandling::SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::vector<double> arrayToWrite;
    // Extract the first element from each nested vector
    for (const auto& innerVector : nestedArrayToWrite) {
        if (!innerVector.empty()) {
            arrayToWrite.push_back(innerVector[0]);
        }
    }

    if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
        if (systemData->printFixedLines) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * systemData->stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if (systemData->printFixedSites) {
            /*outputFileName << (iteration * systemData->stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * systemData->stepsize);
            for (int & fixed_out_val : systemData->fixedOutputSites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (systemData->printAllData) {
        for (int i = 0; i <= systemData->systemTotalSpins; i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * systemData->stepsize) << ","; // Print current time
            else if (i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (systemData->printFixedLines) {
        // iteration >= static_cast<int>(systemData->iterationEnd / 2.0) &&
        if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
            //if (iteration == systemData->iterationEnd) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * systemData->stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;
        }
    } else {
        if (systemData->printAllData) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * systemData->stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
                if (systemData->printFixedSites) {

                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[systemData->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[systemData->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[systemData->systemTotalSpins] << std::endl;

                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[systemData->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[systemData->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(systemData->systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(systemData->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[systemData->systemTotalSpins] << std::endl;
                }
            }
        }
    } */
}

void NMDataHandling::CreateMetadata(bool print_end_time) {

    std::string file_name = "simulation_metadata.txt";

    if (print_end_time) {
        std::ofstream metadata_end;
        metadata_end.open(GV.GetFilePath() + file_name, std::ios_base::app); // append instead of overwrite
        auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        metadata_end << "Finished at:\t" << std::put_time(localtime(&end), "%F %H-%M-%S") << std::endl;
        metadata_end.close();
    }
    else {
        std::ofstream metadata_start(GV.GetFilePath() + file_name);
        CreateFileHeader(metadata_start, "NM 2", true);
        auto start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        metadata_start << "Started at:\t" << std::put_time(localtime(&start), "%F %H-%M-%S") << std::endl;
        metadata_start.close();
    }
}

void NMDataHandling::CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata, int layer) {
    /**
     * Write all non-data information to the output file.
     */
    if (is_metadata) {
        outputFileName << "Key Data\n\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using magDynamics: [" << systemData->useLLG << "]\t\t\t\tUsing Shockwave: [" << systemData->hasShockwave << "]\t\tDrive from LHS: [" << systemData->lhsDrive <<
                       "]\nNumerical Method Used: [" << methodUsed << "]\t\tHas Static Drive: [" << systemData->hasStaticDrive << "]\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0): " << GV.GetStaticBiasField() << " T\t\t\t" << "Dynamic Bias Field (H_D1): " << systemData->dynamicBiasField << " T\n" <<
                          "Dynamic Bias Field Scale Factor: " << systemData->shockwaveInitialStrength << "\t\t" << "Second Dynamic Bias Field (H_D2): " << systemData->shockwaveMax << " T\n" <<
                          "Driving Frequency (f): " << systemData->drivingFreq << "Hz\t\t""Driving Region Start Site: " << systemData->drivingRegionLhs - systemData->numSpinsDamped << "\n" <<
                          "Driving Region End Site: " << systemData->drivingRegionRhs - systemData->numSpinsDamped << " \t\t\t" << "Driving Region Width: " << systemData->drivingRegionWidth << " \n" <<
                          "Max. Sim. Time: " << systemData->maxSimTime << " s\t\t\t\t" << "Min. Exchange Val (J): " << systemData->exchangeEnergyMin  << " T\n" <<
                          "Max. Exchange Val (J): " << systemData->exchangeEnergyMax << " T\t\t\t" << "Max. Iterations: " << systemData->iterationEnd << "\n" <<
                          "No. DataPoints: " << systemData->numberOfDataPoints << " \t\t\t\t" << "No. Spins in Chain: " << systemData->layerSpinsInChain[layer] << "\n" <<
                          "No. Damped Spins: " << systemData->numSpinsDamped << "per side\t\t\t" << "No. Total Spins: " << systemData->layerTotalSpins[layer] << " \n" <<
                          "systemData->stepsize (h): " << systemData->stepsize << "\t\t\t\t" << "Gilbert Damping Factor: " << systemData->gilbertDamping << "\n" <<
                          "Gyromagnetic Ratio (2Pi*Y): " << systemData->gyroMagConst << "\t\t""Shockwave Gradient Time: " << systemData->iterStartShock << "s\n" <<
                          "Shockwave Application Time: " << systemData->shockwaveGradientTime * systemData->stepsize << "s\n" <<
                          std::endl;

        return;
    }
    else {

        outputFileName << "Key Data\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using magDynamics," << systemData->useLLG << ",Using Shockwave," << systemData->hasShockwave << ",Drive from LHS," << systemData->lhsDrive <<
                       ",Numerical Method Used," << methodUsed << ",Has Static Drive," << systemData->hasStaticDrive << "\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                          "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site, Driving Region Width,"
                          "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                          "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, systemData->stepsize (h),Gilbert Damping Factor, Gyromagnetic Ratio (2Pi*Y),"
                          "Shockwave Gradient Time [s], Shockwave Application Time [s]"
                          "\n";

        outputFileName << GV.GetStaticBiasField() << ", " << systemData->dynamicBiasField << ", " << systemData->shockwaveInitialStrength << ", " << systemData->shockwaveMax << ", "
                       << systemData->drivingFreq << ", " << systemData->drivingRegionLhs - systemData->numSpinsDamped << ", " << systemData->drivingRegionRhs - systemData->numSpinsDamped << ", " << systemData->drivingRegionWidth << ", "
                       << systemData->maxSimTime << ", " << systemData->exchangeEnergyMin << ", " << systemData->exchangeEnergyMax << ", " << systemData->iterationEnd << ", " << systemData->numberOfDataPoints << ", "
                       << systemData->layerSpinsInChain[layer] << ", " << systemData->numSpinsDamped << ", " << systemData->layerTotalSpins[layer] << ", " << systemData->stepsize << ", " << systemData->gilbertDamping << ", " << systemData->gyroMagConst << ", "
                       << systemData->iterStartShock << ", " << systemData->shockwaveGradientTime * systemData->stepsize
                       << "\n";

        outputFileName << "\n";
    }

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments );
    outputFileName << "Note(s):," << notesComments << "\n"; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n";

    outputFileName << "\n";

    CreateColumnHeaders(outputFileName, layer);

    std::cout << "\n";
}
void NMDataHandling::CreateColumnHeaders(std::ofstream &outputFileName, int& layer) {
    /**
     * Creates the column headers for each spin site simulated. This code can change often, so compartmentalising it in
     * a separate function is necessary to reduce bugs.
     */
    if (systemData->printAllData or systemData->printFixedLines) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for (int i = 1; i <= systemData->layerTotalSpins[layer]; i++) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if (systemData->printFixedSites) {

        outputFileName << "Time";
        for (int & fixed_out_val : systemData->fixedOutputSites)
            outputFileName << "," << fixed_out_val;
        outputFileName << std::endl;

        //outputFileName << "Time" << ", "
        //               << static_cast<int>(14000) << ","
        //               << static_cast<int>(16000) << ","
        //               << static_cast<int>(18000) << ","
        //               << static_cast<int>(20000) << std::endl;

    }
}
std::vector<double> NMDataHandling::flattenNestedVector(const std::vector<std::vector<double>>& nestedVector) {
    std::vector<double> flattenedVector;

    for (const auto& innerVector : nestedVector) {
        flattenedVector.insert(flattenedVector.end(), innerVector.begin(), innerVector.end());
    }

    return flattenedVector;
}
void NMDataHandling::SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration, int layer) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::vector<double> arrayToWrite;
    // Extract the first element from each nested vector
    for (const auto& innerVector : nestedArrayToWrite) {
        if (!innerVector.empty()) {
            arrayToWrite.push_back(innerVector[0]);
        }
    }

    if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
        if (systemData->printFixedLines) {
            for (int i = 0; i <= systemData->layerTotalSpins[layer]; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * systemData->stepsize) << ",";

                else if (i == systemData->layerTotalSpins[layer])
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if (systemData->printFixedSites) {
            /*outputFileName << (iteration * systemData->stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * systemData->stepsize);
            for (int & fixed_out_val : systemData->fixedOutputSites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (systemData->printAllData) {
        for (int i = 0; i <= systemData->layerTotalSpins[layer]; i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * systemData->stepsize) << ","; // Print current time
            else if (i == systemData->layerTotalSpins[layer])
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (systemData->printFixedLines) {
        // iteration >= static_cast<int>(systemData->iterationEnd / 2.0) &&
        if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
            //if (iteration == systemData->iterationEnd) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * systemData->stepsize) << ",";

                else if (i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;
        }
    } else {
        if (systemData->printAllData) {
            for (int i = 0; i <= systemData->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * systemData->stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (systemData->iterationEnd / systemData->numberOfDataPoints) == 0) {
                if (systemData->printFixedSites) {

                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[systemData->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[systemData->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[systemData->systemTotalSpins] << std::endl;

                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * systemData->stepsize) << ","
                                   << arrayToWrite[systemData->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[systemData->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(systemData->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(systemData->systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(systemData->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[systemData->systemTotalSpins] << std::endl;
                }
            }
        }
    } */
}