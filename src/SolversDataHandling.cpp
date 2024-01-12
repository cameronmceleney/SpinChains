//
// Created by Cameron McEleney on 31/10/2023.
//

// C++ Standard Library
#include "../include/SolversDataHandling.h"

SolversDataHandling::SolversDataHandling( std::shared_ptr<SimulationParameters> sharedSimParams,
                                          std::shared_ptr<SimulationStates> sharedSimStates,
                                          std::shared_ptr<SimulationFlags> sharedSimFlags )

        : SolversSuperClass(std::move(sharedSimParams), std::move(sharedSimStates), std::move(sharedSimFlags)) {}

void SolversDataHandling::CreateFileHeader( std::ofstream &outputFileName, std::string methodUsed, bool is_metadata ) {

    if ( is_metadata ) {
        outputFileName << "Key Data\n\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using LLG: [" << simFlags->shouldUseLLG << "]\t\t\t\tUsing Shockwave: ["
                       << simFlags->hasShockwave << "]\t\tDrive from LHS: [" << simFlags->shouldDriveLHS <<
                       "]\nNumerical Method Used: [" << methodUsed << "]\t\tHas Static Drive: ["
                       << simFlags->isOscillatingZeemanStatic << "]\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0): " << GV.GetStaticBiasField() << " T\t\t\t"
                       << "Dynamic Bias Field (H_D1): " << simParams->oscillatingZeemanStrength << " T\n" <<
                       "Dynamic Bias Field Scale Factor: " << simParams->shockwaveInitialStrength << "\t\t"
                       << "Second Dynamic Bias Field (H_D2): " << simParams->shockwaveMax << " T\n" <<
                       "Driving Frequency (f): " << simParams->drivingFreq << "Hz\t\t""Driving Region Start Site: "
                       << simParams->drivingRegionLhs - simParams->numSpinsInABC << "\n" <<
                       "Driving Region End Site: " << simParams->drivingRegionRhs - simParams->numSpinsInABC
                       << " \t\t\t" << "Driving Region Width: " << simParams->drivingRegionWidth << " \n" <<
                       "Max. Sim. Time: " << simParams->maxSimTime << " s\t\t\t\t" << "Min. Exchange Val (J): "
                       << simParams->exchangeEnergyMin << " T\n" <<
                       "Max. Exchange Val (J): " << simParams->exchangeEnergyMax << " T\t\t\t" << "Max. Iterations: "
                       << simParams->iterationEnd << "\n" <<
                       "No. DataPoints: " << simParams->numberOfDataPoints << " \t\t\t\t" << "No. Spins in Chain: "
                       << simParams->numSpinsInChain << "\n" <<
                       "No. Damped Spins: " << simParams->numSpinsInABC << "per side\t\t\t" << "No. Total Spins: "
                       << simParams->systemTotalSpins << " \n" <<
                       "simParams->stepsize (h): " << simParams->stepsize << "\t\t\t\t" << "Gilbert Damping Factor: "
                       << simParams->gilbertDamping << "\n" <<
                       "Gyromagnetic Ratio (2Pi*Y): " << simParams->gyroMagConst << "\t\t""Shockwave Gradient Time: "
                       << simParams->iterStartShock << "s\n" <<
                       "Shockwave Application Time: " << simParams->shockwaveGradientTime * simParams->stepsize << "s\n"
                       <<
                       std::endl;

        return;
    } else {

        outputFileName << "Key Data\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using magDynamics," << simFlags->shouldUseLLG << ",Using Shockwave,"
                       << simFlags->hasShockwave << ",Drive from LHS," << simFlags->shouldDriveLHS
                       << ",Numerical Method Used," << methodUsed << ",Has Static Drive," << simFlags->isOscillatingZeemanStatic
                       << ",Has Dipolar," << simFlags->hasDipolar << ",Has DMI," << simFlags->hasDMI
                       << ",Has STT," << simFlags->hasSTT << ",Has Zeeman," << simFlags->hasStaticZeeman
                       << ",Has Demag Intense," << simFlags->hasDemagIntense << ",Has Demag FFT," << simFlags->hasDemagFFT
                       << ",Has Shape Anisotropy," << simFlags->hasShapeAnisotropy << "\n";

        outputFileName << "\n";

        outputFileName
                << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                   "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site,Driving Region Width,"
                   "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                   "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, Stepsize (h),Gilbert Damping Factor,Gyromagnetic Ratio (2Pi*Y),"
                   "Shockwave Gradient Time [s],Shockwave Application Time [s],ABC Damping (lower),ABC Damping (upper),"
                   "DMI Constant [T],Saturation Magnetisation [kA/m],Exchange Stiffness [J/m],Anisotropy (Shape) Field [T]"
                   "\n";

        outputFileName << simParams->staticZeemanStrength << ", " << simParams->oscillatingZeemanStrength << ", "
                       << simParams->shockwaveInitialStrength << ", " << simParams->shockwaveMax << ", "
                       << simParams->drivingFreq << ", " << simParams->drivingRegionLhs - simParams->numSpinsInABC << ", "
                       << simParams->drivingRegionRhs - simParams->numSpinsInABC << ", "
                       << simParams->drivingRegionWidth << ", " << simParams->maxSimTime << ", "
                       << simParams->exchangeEnergyMin << ", " << simParams->exchangeEnergyMax << ", "
                       << simParams->iterationEnd << ", " << simParams->numberOfDataPoints << ", "
                       << simParams->numSpinsInChain << ", " << simParams->numSpinsInABC << ", "
                       << simParams->systemTotalSpins << ", " << simParams->stepsize << ", "
                       << simParams->gilbertDamping << ", " << simParams->gyroMagConst << ", "
                       << simParams->iterStartShock << ", " << simParams->shockwaveGradientTime * simParams->stepsize << ", "
                       << simParams->gilbertABCInner << ", " << simParams->gilbertABCOuter << ", "
                       << simParams->dmiConstant << ", " << simParams->satMag << ", "
                       << simParams->exchangeStiffness << ", " << simParams->anisotropyField
                       << "\n";

        outputFileName << "\n";
    }

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments);
    outputFileName << "Note(s):," << notesComments
                   << "\n"; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n";

    outputFileName << "\n";

    CreateColumnHeaders(outputFileName);

    std::cout << "\n";
}

void SolversDataHandling::CreateColumnHeaders( std::ofstream &outputFileName ) {
    /**
     * Creates the column headers for each spin site simulated. This code can change often, so compartmentalising it in
     * a separate function is necessary to reduce bugs.
     */
    if ( simFlags->shouldPrintAllData or simFlags->shouldPrintDiscreteTimes ) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for ( int i = 1; i <= simParams->systemTotalSpins; i++ ) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if ( simFlags->shouldPrintDiscreteSites ) {

        outputFileName << "Time";
        for ( int &fixed_out_val: simStates->fixedOutputSites )
            outputFileName << "," << fixed_out_val;
        outputFileName << std::endl;

        //outputFileName << "Time" << ", "
        //               << static_cast<int>(14000) << ","
        //               << static_cast<int>(16000) << ","
        //               << static_cast<int>(18000) << ","
        //               << static_cast<int>(20000) << std::endl;

    }
}

void SolversDataHandling::InformUserOfCodeType( const std::string &nameNumericalMethod ) {
    /**
     * Informs the user of the code type they are running, including: solver type; special modules.
     */
    if ( simFlags->shouldUseLLG )
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (LLG) code";
    else
        std::cout << "\nYou are running the " << nameNumericalMethod << " Spinchains (Torque) code";

    if ( simFlags->hasShockwave )
        std::cout << " with shockwave module.\n";
    else
        std::cout << ".\n";

}

void SolversDataHandling::PrintVector( std::vector<double> &vectorToPrint, bool shouldExitAfterPrint ) {

    std::cout << "\n\n";

    int count = 0;
    for ( double i: vectorToPrint ) {
        if ( ++count % 10 == 0 )
            std::cout << std::setw(8) << i << std::endl;
        else
            std::cout << std::setw(8) << i << ", ";
    }
    std::cout << "\n\n";

    if ( shouldExitAfterPrint )
        exit(0);
}

void SolversDataHandling::PrintNestedNestedVector( std::vector<std::vector<std::vector<double>>> nestedNestedVector ) {
    // Print the contents of nestedNestedVector1
    std::cout << "{";
    for ( size_t i = 0; i < nestedNestedVector.size(); ++i ) {
        std::cout << "{";
        for ( size_t j = 0; j < nestedNestedVector[i].size(); ++j ) {
            std::cout << "{";
            for ( size_t k = 0; k < nestedNestedVector[i][j].size(); ++k ) {
                std::cout << nestedNestedVector[i][j][k];
                if ( k < nestedNestedVector[i][j].size() - 1 ) {
                    std::cout << ",";
                }
            }
            std::cout << "}";
            if ( j < nestedNestedVector[i].size() - 1 ) {
                std::cout << ",";
            }
        }
        std::cout << "}";
        if ( i < nestedNestedVector.size() - 1 ) {
            std::cout << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}

void SolversDataHandling::SaveDataToFile( std::ofstream &outputFileName, std::vector<double> &arrayToWrite,
                                          int &iteration ) {
    std::cout.precision(6);
    std::cout << std::scientific;

    if ( iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0 ) {
        if ( simFlags->shouldPrintDiscreteTimes ) {
            for ( int i = 0; i <= simParams->systemTotalSpins; i++ ) {
                // Steps through vectors containing all mag. moment components and saves to files
                if ( i == 0 )
                    // Print current time
                    outputFileName << (iteration * simParams->stepsize) << ",";

                else if ( i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if ( simFlags->shouldPrintDiscreteSites ) {
            /*outputFileName << (iteration * simParams->stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * simParams->stepsize);
            for ( int &fixed_out_val: simStates->fixedOutputSites )
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if ( simFlags->shouldPrintAllData ) {
        for ( int i = 0; i <= simParams->systemTotalSpins; i++ ) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if ( i == 0 )
                outputFileName << (iteration * simParams->stepsize) << ","; // Print current time
            else if ( i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (simFlags->shouldPrintDiscreteTimes) {
        // iteration >= static_cast<int>(simParams->iterationEnd / 2.0) &&
        if (iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0) {
            //if (iteration == simParams->iterationEnd) {
            for (int i = 0; i <= simParams->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * simParams->stepsize) << ",";

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
        if (simFlags->shouldPrintAllData) {
            for (int i = 0; i <= simParams->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * simParams->stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0) {
                if (simFlags->shouldPrintDiscreteSites) {

                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[simParams->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[simParams->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[simParams->systemTotalSpins] << std::endl;

                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[simParams->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[simParams->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(simParams->systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(simParams->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[simParams->systemTotalSpins] << std::endl;
                }
            }
        }
    } */
}

void SolversDataHandling::SaveDataToFileMultilayer( std::ofstream &outputFileName,
                                                    std::vector<std::vector<double>> &nestedArrayToWrite,
                                                    int &iteration ) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::vector<double> arrayToWrite;
    // Extract the first element from each nested vector
    for ( const auto &innerVector: nestedArrayToWrite ) {
        if ( !innerVector.empty()) {
            arrayToWrite.push_back(innerVector[0]);
        }
    }

    if ( iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0 ) {
        if ( simFlags->shouldPrintDiscreteTimes ) {
            for ( int i = 0; i <= simParams->systemTotalSpins; i++ ) {
                // Steps through vectors containing all mag. moment components and saves to files
                if ( i == 0 )
                    // Print current time
                    outputFileName << (iteration * simParams->stepsize) << ",";

                else if ( i == GV.GetNumSpins())
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if ( simFlags->shouldPrintDiscreteSites ) {
            /*outputFileName << (iteration * simParams->stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * simParams->stepsize);
            for ( int &fixed_out_val: simStates->fixedOutputSites )
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if ( simFlags->shouldPrintAllData ) {
        for ( int i = 0; i <= simParams->systemTotalSpins; i++ ) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if ( i == 0 )
                outputFileName << (iteration * simParams->stepsize) << ","; // Print current time
            else if ( i == GV.GetNumSpins())
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (simFlags->shouldPrintDiscreteTimes) {
        // iteration >= static_cast<int>(simParams->iterationEnd / 2.0) &&
        if (iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0) {
            //if (iteration == simParams->iterationEnd) {
            for (int i = 0; i <= simParams->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * simParams->stepsize) << ",";

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
        if (simFlags->shouldPrintAllData) {
            for (int i = 0; i <= simParams->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * simParams->stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0) {
                if (simFlags->shouldPrintDiscreteSites) {

                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[simParams->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[simParams->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[simParams->systemTotalSpins] << std::endl;

                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[simParams->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[simParams->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(simParams->systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(simParams->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[simParams->systemTotalSpins] << std::endl;
                }
            }
        }
    } */
}

void SolversDataHandling::CreateMetadata( bool print_end_time ) {

    std::string file_name = "simulation_metadata.txt";

    if ( print_end_time ) {
        std::ofstream metadata_end;
        metadata_end.open(GV.GetFilePath() + file_name, std::ios_base::app); // append instead of overwrite
        auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        metadata_end << "Finished at:\t" << std::put_time(localtime(&end), "%F %H-%M-%S") << std::endl;
        metadata_end.close();
    } else {
        std::ofstream metadata_start(GV.GetFilePath() + file_name);
        CreateFileHeader(metadata_start, "NM 2", true);
        auto start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        metadata_start << "Started at:\t" << std::put_time(localtime(&start), "%F %H-%M-%S") << std::endl;
        metadata_start.close();
    }
}

void SolversDataHandling::CreateFileHeader( std::ofstream &outputFileName, std::string methodUsed, bool is_metadata,
                                            int layer ) {
    /**
     * Write all non-data information to the output file.
     */
    if ( is_metadata ) {
        outputFileName << "Key Data\n\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using magDynamics: [" << simFlags->shouldUseLLG << "]\t\t\t\tUsing Shockwave: ["
                       << simFlags->hasShockwave << "]\t\tDrive from LHS: [" << simFlags->shouldDriveLHS <<
                       "]\nNumerical Method Used: [" << methodUsed << "]\t\tHas Static Drive: ["
                       << simFlags->isOscillatingZeemanStatic << "]\n";

        outputFileName << "\n";

        outputFileName << "Static Bias Field (H0): " << GV.GetStaticBiasField() << " T\t\t\t"
                       << "Dynamic Bias Field (H_D1): " << simParams->oscillatingZeemanStrength << " T\n" <<
                       "Dynamic Bias Field Scale Factor: " << simParams->shockwaveInitialStrength << "\t\t"
                       << "Second Dynamic Bias Field (H_D2): " << simParams->shockwaveMax << " T\n" <<
                       "Driving Frequency (f): " << simParams->drivingFreq << "Hz\t\t""Driving Region Start Site: "
                       << simParams->drivingRegionLhs - simParams->numSpinsInABC << "\n" <<
                       "Driving Region End Site: " << simParams->drivingRegionRhs - simParams->numSpinsInABC
                       << " \t\t\t" << "Driving Region Width: " << simParams->drivingRegionWidth << " \n" <<
                       "Max. Sim. Time: " << simParams->maxSimTime << " s\t\t\t\t" << "Min. Exchange Val (J): "
                       << simParams->exchangeEnergyMin << " T\n" <<
                       "Max. Exchange Val (J): " << simParams->exchangeEnergyMax << " T\t\t\t" << "Max. Iterations: "
                       << simParams->iterationEnd << "\n" <<
                       "No. DataPoints: " << simParams->numberOfDataPoints << " \t\t\t\t" << "No. Spins in Chain: "
                       << simStates->layerSpinsInChain[layer] << "\n" <<
                       "No. Damped Spins: " << simParams->numSpinsInABC << "per side\t\t\t" << "No. Total Spins: "
                       << simStates->layerTotalSpins[layer] << " \n" <<
                       "Stepsize (h): " << simParams->stepsize << "\t\t\t\t" << "Gilbert Damping Factor (Chain: "
                       << simParams->gilbertDamping << "\n" <<
                       "Gyromagnetic Ratio (2Pi*Y): " << simParams->gyroMagConst << "\t\t""Shockwave Gradient Time: "
                       << simParams->iterStartShock << "s\n" <<
                       "Shockwave Application Time: " << simParams->shockwaveGradientTime * simParams->stepsize << "s\n"
                       <<
                       std::endl;

        return;
    } else {

        outputFileName << "Key Data\n";

        outputFileName << "[Booleans where (1) indicates (True) and (0) indicates (False)]\n";

        outputFileName << "Using magDynamics," << simFlags->shouldUseLLG << ",Using Shockwave,"
                       << simFlags->hasShockwave << ",Drive from LHS," << simFlags->shouldDriveLHS <<
                       ",Numerical Method Used," << methodUsed << ",Has Static Drive," << simFlags->isOscillatingZeemanStatic
                       << "\n";

        outputFileName << "\n";

        outputFileName
                << "Static Bias Field (H0) [T],Dynamic Bias Field (H_D1) [T],Dynamic Bias Field Scale Factor,Second Dynamic Bias Field (H_D2)[T],"
                   "Driving Frequency (f) [Hz],Driving Region Start Site,Driving Region End Site, Driving Region Width,"
                   "Max. Sim. Time [s],Min. Exchange Val (J)[T],Max. Exchange Val (J)[T],Max. Iterations,No. DataPoints,"
                   "No. Spins in Chain (N),No. Damped Spins (per side),No. Total Spins, simParams->stepsize (h),Gilbert Damping Factor, Gyromagnetic Ratio (2Pi*Y),"
                   "Shockwave Gradient Time [s], Shockwave Application Time [s]"
                   "\n";

        outputFileName << GV.GetStaticBiasField() << ", " << simParams->oscillatingZeemanStrength << ", "
                       << simParams->shockwaveInitialStrength << ", " << simParams->shockwaveMax << ", "
                       << simParams->drivingFreq << ", " << simParams->drivingRegionLhs - simParams->numSpinsInABC
                       << ", " << simParams->drivingRegionRhs - simParams->numSpinsInABC << ", "
                       << simParams->drivingRegionWidth << ", "
                       << simParams->maxSimTime << ", " << simParams->exchangeEnergyMin << ", "
                       << simParams->exchangeEnergyMax << ", " << simParams->iterationEnd << ", "
                       << simParams->numberOfDataPoints << ", "
                       << simStates->layerSpinsInChain[layer] << ", " << simParams->numSpinsInABC << ", "
                       << simStates->layerTotalSpins[layer] << ", " << simParams->stepsize << ", "
                       << simParams->gilbertDamping << ", " << simParams->gyroMagConst << ", "
                       << simParams->iterStartShock << ", " << simParams->shockwaveGradientTime * simParams->stepsize
                       << "\n";

        outputFileName << "\n";
    }

    std::string notesComments;
    std::cout << "Enter any notes for this simulation: ";
    std::cin.ignore();
    std::getline(std::cin, notesComments);
    outputFileName << "Note(s):," << notesComments
                   << "\n"; // Adding comma ensures the note itself is in a different csv cell to the term 'Note(s):'

    outputFileName << "[Column heading indicates the spin site (#) being recorded. Data is for the (mx) component]\n";

    outputFileName << "\n";

    CreateColumnHeaders(outputFileName, layer);

    std::cout << "\n";
}

void SolversDataHandling::CreateColumnHeaders( std::ofstream &outputFileName, int &layer ) {
    /**
     * Creates the column headers for each spin site simulated. This code can change often, so compartmentalising it in
     * a separate function is necessary to reduce bugs.
     */
    if ( simFlags->shouldPrintAllData or simFlags->shouldPrintDiscreteTimes ) {
        // Print column heading for every spin simulated.
        outputFileName << "Time [s], ";
        for ( int i = 1; i <= simStates->layerTotalSpins[layer]; i++ ) {
            outputFileName << i << ", ";
        }
        outputFileName << std::endl;

    } else if ( simFlags->shouldPrintDiscreteSites ) {

        outputFileName << "Time";
        for ( int &fixed_out_val: simStates->fixedOutputSites )
            outputFileName << "," << fixed_out_val;
        outputFileName << std::endl;

        //outputFileName << "Time" << ", "
        //               << static_cast<int>(14000) << ","
        //               << static_cast<int>(16000) << ","
        //               << static_cast<int>(18000) << ","
        //               << static_cast<int>(20000) << std::endl;

    }
}

std::vector<double> SolversDataHandling::flattenNestedVector( const std::vector<std::vector<double>> &nestedVector ) {
    std::vector<double> flattenedVector;

    for ( const auto &innerVector: nestedVector ) {
        flattenedVector.insert(flattenedVector.end(), innerVector.begin(), innerVector.end());
    }

    return flattenedVector;
}

void SolversDataHandling::SaveDataToFileMultilayer( std::ofstream &outputFileName,
                                                    std::vector<std::vector<double>> &nestedArrayToWrite,
                                                    int &iteration, int layer ) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::vector<double> arrayToWrite;
    // Extract the first element from each nested vector
    for ( const auto &innerVector: nestedArrayToWrite ) {
        if ( !innerVector.empty()) {
            arrayToWrite.push_back(innerVector[0]);
        }
    }

    if ( iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0 ) {
        if ( simFlags->shouldPrintDiscreteTimes ) {
            for ( int i = 0; i <= simStates->layerTotalSpins[layer]; i++ ) {
                // Steps through vectors containing all mag. moment components and saves to files
                if ( i == 0 )
                    // Print current time
                    outputFileName << (iteration * simParams->stepsize) << ",";

                else if ( i == simStates->layerTotalSpins[layer] )
                    // Ensures that the final line doesn't contain a comma.
                    outputFileName << arrayToWrite[i] << std::flush;

                else
                    // For non-special values, write the data.
                    outputFileName << arrayToWrite[i] << ", ";
            }
            // Take new line after current row is finished being written.
            outputFileName << std::endl;

            return;
        } else if ( simFlags->shouldPrintDiscreteSites ) {
            /*outputFileName << (iteration * simParams->stepsize) << ","
               << arrayToWrite[14000] << ","
               << arrayToWrite[16000] << ","
               << arrayToWrite[18000] << ","
               << arrayToWrite[20000] << std::endl;
               */
            outputFileName << (iteration * simParams->stepsize);
            for ( int &fixed_out_val: simStates->fixedOutputSites )
                outputFileName << "," << arrayToWrite[fixed_out_val];
            outputFileName << std::endl;

            return;
        }
    }

    if ( simFlags->shouldPrintAllData ) {
        for ( int i = 0; i <= simStates->layerTotalSpins[layer]; i++ ) {
            // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
            if ( i == 0 )
                outputFileName << (iteration * simParams->stepsize) << ","; // Print current time
            else if ( i == simStates->layerTotalSpins[layer] )
                outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
            else
                outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
        }
        outputFileName << std::endl;

        return;
    }

    /*
    if (simFlags->shouldPrintDiscreteTimes) {
        // iteration >= static_cast<int>(simParams->iterationEnd / 2.0) &&
        if (iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0) {
            //if (iteration == simParams->iterationEnd) {
            for (int i = 0; i <= simParams->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components and saves to files
                if (i == 0)
                    // Print current time
                    outputFileName << (iteration * simParams->stepsize) << ",";

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
        if (simFlags->shouldPrintAllData) {
            for (int i = 0; i <= simParams->systemTotalSpins; i++) {
                // Steps through vectors containing all mag. moment components found at the end of RK2-Stage 2, and saves to files
                if (i == 0)
                    outputFileName << (iteration * simParams->stepsize) << ","; // Print current time
                else if (i == GV.GetNumSpins())
                    outputFileName << arrayToWrite[i] << std::flush; // Ensures that the final line doesn't contain a comma.
                else
                    outputFileName << arrayToWrite[i] << ","; // For non-special values, write the data.
            }
            outputFileName << std::endl; // Take new line after current row is finished being written.
        } else {
            if (iteration % (simParams->iterationEnd / simParams->numberOfDataPoints) == 0) {
                if (simFlags->shouldPrintDiscreteSites) {

                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[simParams->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[simParams->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(1500)] << ","
                                   << arrayToWrite[static_cast<int>(2500)] << ","
                                   << arrayToWrite[static_cast<int>(3500)] << ","
                                   << arrayToWrite[simParams->systemTotalSpins] << std::endl;

                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[400] << ","
                                   << arrayToWrite[1500] << ","
                                   << arrayToWrite[3000] << ","
                                   << arrayToWrite[4500] << ","
                                   << arrayToWrite[5600] << std::endl;
                } else {
                    outputFileName << (iteration * simParams->stepsize) << ","
                                   << arrayToWrite[simParams->drivingRegionLhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->drivingRegionWidth / 2.0)] << ","
                                   << arrayToWrite[simParams->drivingRegionRhs] << ","
                                   << arrayToWrite[static_cast<int>(simParams->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[static_cast<int>(simParams->systemTotalSpins / 2.0)] << ","
                                   << arrayToWrite[3 * static_cast<int>(simParams->systemTotalSpins / 4.0)] << ","
                                   << arrayToWrite[simParams->systemTotalSpins] << std::endl;
                }
            }
        }
    } */
}