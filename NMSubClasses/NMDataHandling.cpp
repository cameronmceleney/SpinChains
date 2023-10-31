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

        outputFileName << "Using LLG: [" << useLLG << "]\t\t\t\tUsing Shockwave: [" << hasShockwave << "]\t\tDrive from LHS: [" << lhsDrive <<
                       "]\nNumerical Method Used: [" << methodUsed << "]\t\tHas Static Drive: [" << hasStaticDrive << "]\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0): " << GV.GetStaticBiasField() << " T\t\t\t" << "Dynamic Bias Field (H_D1): " << dynamicBiasField << " T\n" <<
                          "Dynamic Bias Field Scale Factor: " << shockwaveInitialStrength << "\t\t" << "Second Dynamic Bias Field (H_D2): " << shockwaveMax << " T\n" <<
                          "Driving Frequency (f): " << drivingFreq << "Hz\t\t""Driving Region Start Site: " << drivingRegionLhs - numSpinsDamped << "\n" <<
                          "Driving Region End Site: " << drivingRegionRhs - numSpinsDamped << " \t\t\t" << "Driving Region Width: " << drivingRegionWidth << " \n" <<
                          "Max. Sim. Time: " << maxSimTime << " s\t\t\t\t" << "Min. Exchange Val (J): " << exchangeEnergyMin  << " T\n" <<
                          "Max. Exchange Val (J): " << exchangeEnergyMax << " T\t\t\t" << "Max. Iterations: " << iterationEnd << "\n" <<
                          "No. DataPoints: " << numberOfDataPoints << " \t\t\t\t" << "No. Spins in Chain: " << numSpinsInChain << "\n" <<
                          "No. Damped Spins: " << numSpinsDamped << "per side\t\t\t" << "No. Total Spins: " << systemTotalSpins << " \n" <<
                          "Stepsize (h): " << stepsize << "\t\t\t\t" << "Gilbert Damping Factor: " << gilbertDamping << "\n" <<
                          "Gyromagnetic Ratio (2Pi*Y): " << gyroMagConst << "\t\t""Shockwave Gradient Time: " << iterStartShock << "s\n" <<
                          "Shockwave Application Time: " << shockwaveGradientTime * stepsize << "s\n" <<
                          std::endl;

        return;
    }
    else {

        outputFileName << "Key Data\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using LLG," << useLLG << ",Using Shockwave," << hasShockwave << ",Drive from LHS," << lhsDrive <<
                       ",Numerical Method Used," << methodUsed << ",Has Static Drive," << hasStaticDrive << "\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                          "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site, Driving Region Width,"
                          "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                          "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, Stepsize (h),Gilbert Damping Factor, Gyromagnetic Ratio (2Pi*Y),"
                          "Shockwave Gradient Time [s], Shockwave Application Time [s]"
                          "\n";

        outputFileName << GV.GetStaticBiasField() << ", " << dynamicBiasField << ", " << shockwaveInitialStrength << ", " << shockwaveMax << ", "
                       << drivingFreq << ", " << drivingRegionLhs - numSpinsDamped << ", " << drivingRegionRhs - numSpinsDamped << ", " << drivingRegionWidth << ", "
                       << maxSimTime << ", " << exchangeEnergyMin << ", " << exchangeEnergyMax << ", " << iterationEnd << ", " << numberOfDataPoints << ", "
                       << numSpinsInChain << ", " << numSpinsDamped << ", " << systemTotalSpins << ", " << stepsize << ", " << gilbertDamping << ", " << gyroMagConst << ", "
                       << iterStartShock << ", " << shockwaveGradientTime * stepsize
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
    if (printAllData or printFixedLines) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for (int i = 1; i <= systemTotalSpins; i++) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if (printFixedSites) {

        outputFileName << "Time";
        for (int & fixed_out_val : fixedOutputSites)
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
    if (useLLG)
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (LLG) code";
    else
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (Torque) code";

    if (hasShockwave)
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

    if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
        if (printFixedLines) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * stepsize) << ",";

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
        } else if (printFixedSites) {
            /*outputFileName << (iteration * stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * stepsize);
            for (int & fixed_out_val : fixedOutputSites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (printAllData) {
        for (int i = 0; i <= systemTotalSpins; i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * stepsize) << ","; // Print current time
            else if (i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (printFixedLines) {
        // iteration >= static_cast<int>(iterationEnd / 2.0) &&
        if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
            //if (iteration == iterationEnd) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * stepsize) << ",";

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
        if (printAllData) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
                if (printFixedSites) {

                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[systemTotalSpins] << std::endl;

                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[systemTotalSpins] << std::endl;
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

    if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
        if (printFixedLines) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * stepsize) << ",";

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
        } else if (printFixedSites) {
            /*outputFileName << (iteration * stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * stepsize);
            for (int & fixed_out_val : fixedOutputSites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (printAllData) {
        for (int i = 0; i <= systemTotalSpins; i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * stepsize) << ","; // Print current time
            else if (i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (printFixedLines) {
        // iteration >= static_cast<int>(iterationEnd / 2.0) &&
        if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
            //if (iteration == iterationEnd) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * stepsize) << ",";

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
        if (printAllData) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
                if (printFixedSites) {

                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[systemTotalSpins] << std::endl;

                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[systemTotalSpins] << std::endl;
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

        outputFileName << "Using LLG: [" << useLLG << "]\t\t\t\tUsing Shockwave: [" << hasShockwave << "]\t\tDrive from LHS: [" << lhsDrive <<
                       "]\nNumerical Method Used: [" << methodUsed << "]\t\tHas Static Drive: [" << hasStaticDrive << "]\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0): " << GV.GetStaticBiasField() << " T\t\t\t" << "Dynamic Bias Field (H_D1): " << dynamicBiasField << " T\n" <<
                          "Dynamic Bias Field Scale Factor: " << shockwaveInitialStrength << "\t\t" << "Second Dynamic Bias Field (H_D2): " << shockwaveMax << " T\n" <<
                          "Driving Frequency (f): " << drivingFreq << "Hz\t\t""Driving Region Start Site: " << drivingRegionLhs - numSpinsDamped << "\n" <<
                          "Driving Region End Site: " << drivingRegionRhs - numSpinsDamped << " \t\t\t" << "Driving Region Width: " << drivingRegionWidth << " \n" <<
                          "Max. Sim. Time: " << maxSimTime << " s\t\t\t\t" << "Min. Exchange Val (J): " << exchangeEnergyMin  << " T\n" <<
                          "Max. Exchange Val (J): " << exchangeEnergyMax << " T\t\t\t" << "Max. Iterations: " << iterationEnd << "\n" <<
                          "No. DataPoints: " << numberOfDataPoints << " \t\t\t\t" << "No. Spins in Chain: " << layerSpinsInChain[layer] << "\n" <<
                          "No. Damped Spins: " << numSpinsDamped << "per side\t\t\t" << "No. Total Spins: " << layerTotalSpins[layer] << " \n" <<
                          "Stepsize (h): " << stepsize << "\t\t\t\t" << "Gilbert Damping Factor: " << gilbertDamping << "\n" <<
                          "Gyromagnetic Ratio (2Pi*Y): " << gyroMagConst << "\t\t""Shockwave Gradient Time: " << iterStartShock << "s\n" <<
                          "Shockwave Application Time: " << shockwaveGradientTime * stepsize << "s\n" <<
                          std::endl;

        return;
    }
    else {

        outputFileName << "Key Data\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using LLG," << useLLG << ",Using Shockwave," << hasShockwave << ",Drive from LHS," << lhsDrive <<
                       ",Numerical Method Used," << methodUsed << ",Has Static Drive," << hasStaticDrive << "\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                          "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site, Driving Region Width,"
                          "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                          "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, Stepsize (h),Gilbert Damping Factor, Gyromagnetic Ratio (2Pi*Y),"
                          "Shockwave Gradient Time [s], Shockwave Application Time [s]"
                          "\n";

        outputFileName << GV.GetStaticBiasField() << ", " << dynamicBiasField << ", " << shockwaveInitialStrength << ", " << shockwaveMax << ", "
                       << drivingFreq << ", " << drivingRegionLhs - numSpinsDamped << ", " << drivingRegionRhs - numSpinsDamped << ", " << drivingRegionWidth << ", "
                       << maxSimTime << ", " << exchangeEnergyMin << ", " << exchangeEnergyMax << ", " << iterationEnd << ", " << numberOfDataPoints << ", "
                       << layerSpinsInChain[layer] << ", " << numSpinsDamped << ", " << layerTotalSpins[layer] << ", " << stepsize << ", " << gilbertDamping << ", " << gyroMagConst << ", "
                       << iterStartShock << ", " << shockwaveGradientTime * stepsize
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
    if (printAllData or printFixedLines) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for (int i = 1; i <= layerTotalSpins[layer]; i++) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if (printFixedSites) {

        outputFileName << "Time";
        for (int & fixed_out_val : fixedOutputSites)
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

    if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
        if (printFixedLines) {
            for (int i = 0; i <= layerTotalSpins[layer]; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * stepsize) << ",";

                else if (i == layerTotalSpins[layer])
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if (printFixedSites) {
            /*outputFileName << (iteration * stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * stepsize);
            for (int & fixed_out_val : fixedOutputSites)
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if (printAllData) {
        for (int i = 0; i <= layerTotalSpins[layer]; i++) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if (i == 0)
                outputFileName << (iteration * stepsize) << ","; // Print current time
            else if (i == layerTotalSpins[layer])
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (printFixedLines) {
        // iteration >= static_cast<int>(iterationEnd / 2.0) &&
        if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
            //if (iteration == iterationEnd) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * stepsize) << ",";

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
        if (printAllData) {
            for (int i = 0; i <= systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (iterationEnd / numberOfDataPoints) == 0) {
                if (printFixedSites) {

                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[systemTotalSpins] << std::endl;

                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * stepsize) << ","
                                   << arrayToWrite[drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[systemTotalSpins] << std::endl;
                }
            }
        }
    } */
}