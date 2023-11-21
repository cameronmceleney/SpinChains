//
// Created by Cameron Aidan McEleney on 20/11/2023.
//

#ifndef SPINCHAINS_SIMULATIONFLAGS_H
#define SPINCHAINS_SIMULATIONFLAGS_H

class SimulationFlags {
public:
    bool                centralDrive;                               // Drive from the centre of the chain if (true)
    bool                driveAllLayers;
    bool                dualDrive;                                  // Drive from both sides of the system
    bool                hasShockwave;                               // Simulation contains a single driving bias field if (false).

    bool                hasStaticDrive;                             // Selects (if true) whether drive has sinusoidal term
    bool                isShockwaveOn;                              // Tests if the conditions to trigger a shockwave have been reached. Not to be altered by the user.
    bool                isShockwaveAtMax;                           // Tests if the shockwave is at its maximum amplitude. Not to be altered by the user.

    bool                lhsDrive;                                   // Drive from the LHS
    bool                rhsDrive;                                   // Drive from the RHS
    bool                printAllData;                               // Saves the m-component(s) of every spin at every iteration. WARNING: leads to huge output files.
    bool                printFixedLines;                            // Saves m-component(s) of every spin at regular intervals. Total save points are set by numberOfDataPoints.
    bool                printFixedSites;                            // Saves a discrete set of m-component(s) at regular intervals governed by numberOfDataPoints.

    bool                shouldDriveCease;                           // Internal flag to indicate if the driving field should cut off at a given time.
    bool                shouldTrackMValues;                         // Monitor the norm of all the m-values; if approx. 1.0 then the error is likely to be massive; discard that dataset.
    bool                useLLG;                                     // Uses the Torque equation components if (false).
    bool                useSLLG;

    bool                useDipolar;
    bool                useZeeman;
    bool                useDemagIntense;                            // todo doesn't work
    bool                useDemagFft;                                // todo doesn't work
    bool                useMultilayer;
    bool                debugFunc;

    bool isFm;
};

#endif //SPINCHAINS_SIMULATIONFLAGS_H
