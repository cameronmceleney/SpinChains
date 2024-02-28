//
// Created by Cameron Aidan McEleney on 20/11/2023.
//

#ifndef SPINCHAINS_SIMULATIONFLAGS_H
#define SPINCHAINS_SIMULATIONFLAGS_H

class SimulationFlags {
public:
    bool                shouldDriveCentre;                               // Drive from the centre of the chain if (true)
    bool                shouldDriveAllLayers;
    bool                shouldDriveBothSides;                                  // Drive from both sides of the system
    bool                hasShockwave;                               // Simulation contains a single driving bias field if (false).

    bool                isOscillatingZeemanStatic;                             // Selects (if true) whether drive has sinusoidal term
    bool                isShockwaveOn;                              // Tests if the conditions to trigger a shockwave have been reached. Not to be altered by the user.
    bool                isShockwaveAtMax;                           // Tests if the shockwave is at its maximum amplitude. Not to be altered by the user.

    bool                shouldDriveLHS;                                   // Drive from the LHS
    bool                shouldDriveRHS;                                   // Drive from the RHS
    bool                shouldDriveDiscreteSites;
    bool                hasCustomDrivePosition;
    bool                hasDemag1DThinFilm;

    bool                shouldPrintAllData;                               // Saves the m-component(s) of every spin at every iteration. WARNING: leads to huge output files.
    bool                shouldPrintDiscreteTimes;                            // Saves m-component(s) of every spin at regular intervals. Total save points are set by numberOfDataPoints.
    bool                shouldPrintDiscreteSites;                            // Saves a discrete set of m-component(s) at regular intervals governed by numberOfDataPoints.

    bool                shouldDriveCeaseEarly;                           // Internal flag to indicate if the driving field should cut off at a given time.
    bool                shouldTrackMagneticMomentNorm;                         // Monitor the norm of all the m-values; if approx. 1.0 then the error is likely to be massive; discard that dataset.
    bool                shouldUseLLG;                                     // Uses the Torque equation components if (false).
    bool                shouldUseSLLG;

    bool                hasDipolar;
    bool                hasDMI;
    bool                hasSTT;
    bool                hasStaticZeeman;
    bool                hasDemagIntense;                            // todo doesn't work
    bool                hasDemagFFT;                                // todo doesn't work
    bool                hasMultipleStackedLayers;
    bool                hasMultipleAppendedLayers;
    bool                debugFunc;
    bool                hasShapeAnisotropy;
    int                 preferredDirection;                         // Establish preferred magnetic field direction; set static external field and DMI along this in 1D. [0] = x, [1] = y, [2] = z

    bool                isFerromagnetic;
    bool                resetSimState = false;
    bool                hasSingleExchangeRegion = false;

    bool                hasGradientWithinDrivingRegion;
};

#endif //SPINCHAINS_SIMULATIONFLAGS_H
