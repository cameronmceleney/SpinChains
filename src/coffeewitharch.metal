#include <metal_stdlib>
using namespace metal;

/*
 * Would be good to turn these structs into function constants (to control the code paths to compile)
 *
 * e.g.
 *  struct VertexOutput {
 *      float4 position [[position]];
 *      float4 color;
 *  }
 *
 *  #ifdef OFFSET_DEFINED
 *      vOut.position += vIn.position;
 *  #endif
 *  would become
 *
 *  struct VertexOutput {
 *      float4 position [[position]];
 *      float4 color;
 *  }
 *  constant bool offset_defined [[function_constant(0)]];
 *  if (offset_defined) { vOut.position += vIn.offset;
 *
 */

typedef unsigned short u_short;
typedef unsigned int u_int;

struct SimulationParams {
    // todo separate this out into several compositions (flags, data structures, etc)
    float PERMEABILITY_IRON = 2.22; // cobalt = 1.72, iron = 2.22, nickel = 0.6;

    // Sites to be printed if shouldPrintDiscreteSites is TRUE.

    float               gpu_current_time;
    float               ambientTemperature;
    float               anisotropyField;
    float              staticZeemanStrength;

    float               drivingAngFreq;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
    float               drivingFreq;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
    u_short              drivingRegionLhs;                         // The position of the spin which is leftmost in the driving region.
    u_short                 drivingRegionRhs;                         // The position of the spin which is rightmost in the driving region.

    u_int                 drivingRegionWidth;                       // Driving region width.
    float              oscillatingZeemanStrength;                         // Driving field amplitude [T] (caution: papers often give in [mT]).
    u_int                 forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
    float              dipoleConstant;                           // Scaling factor which is constant across dipolar interaction calculations.
    float              dmiConstant;
    float              exchangeStiffness;

    float              gilbertDamping;                           // Gilbert damping factor for main chain.
    float              gilbertABCInner;                          // The lower Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the chain meets the ABC.
    float              gilbertABCOuter;                          // The upper Gilbert damping factor for the Absorbing boundary conditions (ABCs) at the point where the ABC meets the pinned sites.
    float              gyroMagConst;                             // Gyromagnetic ratio of an electron [GHz/T].
    float              latticeConstant;

    float              spinPolarisation;                         // Spin polarisation of the spin current.
    float              spinTransferEfficiency;                   // Spin transfer efficiency of the spin current.

    u_int                 iterationEnd;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    u_int                 iterationStart = 0;                        // The iteration step that the program will begin at (Default: 0.0)
    float              iterStartShock;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the maxSimTime.

    float              iterEndShock;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the maxSimTime.
    float              largestMNorm = 1e-40;                     // Computes sqrt(_mx**2 + _my**2 + _mz**2) for each site at each moment to track any abnormalities. Initialises to be arbitrarily small
    float              maxSimTime;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [mxInit + myInit + mzInit]  CANNOT sum to greater than 1.0
    u_short                 numNeighbours;
    u_int                 numberOfDataPoints;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    u_int                 numberOfSpinPairs;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    u_int                 numSpinsInABC;                           // Number of spins in the damped regions (previously called _numGilbert).

    u_int                numSpinsInChain;                          // The number of spin sites in the spin chain to be simulated.
    u_int                 systemTotalSpins;                         // The total number of spins in the system (chain plus ABCs).
    float              satMag;                                   // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)

    float              shockwaveGradientTime;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    float              shockwaveInitialStrength;                 // Initial strength of the shockwave before shockwaveScaling occurs. (Default: = oscillatingZeemanStrength)
    float              shockwaveMax;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    float              shockwaveScaling;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving

    float              shockwaveStepsize;                        // Size of incremental increase in shockwave amplitude.
    float              stepsize;                                 // Stepsize between values
    float              stepsizeHalf;                             // Separately defined to avoid repeated unnecessary calculations inside loops
    u_short                 numLayers;
    float              recordingInterval;                        // Time between each data point recorded.
    u_short                 layerOfInterest;

    float exchangeEnergyMin;
    float exchangeEnergyMax;

    float totalTime = 0;

    u_short numSpinsDRPeak, numSpinsDRGradient;
    u_short dmiRegionLhs, dmiRegionRhs, numSpinsDmiPeak, numSpinsDmiGradient, numSpinsDmiWidth, dmiRegionOffset;
    u_short dampingRegionLhs, dampingRegionRhs, numSpinsDampingPeak, numSpinsDampingGradient, numSpinsDampingWidth, dampingRegionOffset;

    float dampingGradientPeak, dampingGradientGradient;

};

struct SimulationFlags {
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
    bool                resetSimState;
    bool                hasSingleExchangeRegion;

    bool                hasGradientRegionForOscillatingZeeman;
    bool                hasGradientRegionForDmi;
    bool                hasGradientRegionForDamping;

    bool                shouldDmiGradientMirrorOscillatingZeeman;
    bool                shouldDampingGradientMirrorOscillatingZeeman;

    bool                shouldRestrictDmiToWithinGradientRegion;
    bool                isOscillatingZeemanLinearAcrossMap;
    bool                isDmiLinearAcrossMap;
    bool                isDampingLinearAcrossMap;

    bool                useGenerateABCUpdated;
    bool                forceSequentialOperation;
};

// Function constants are a cleaner way to achieve these Boolean checks rather than using a struct
// constant bool hasGradientRegionForDmi [[ function_constant(0) ]];
// constant bool hasGradientRegionForOscillatingZeeman [[ function_constant(1) ]];
// constant bool shouldRestrictDmiToWithinGradientRegion [[ function_constant(2) ]];
// constant bool shouldUseLLG [[ function_constant(3) ]];
// constant bool shouldDriveAllLayers [[ function_constant(4) ]];
// constant bool isOscillatingZeemanStatic [[ function_constant(5) ]];
// constant bool shouldDriveDiscreteSites [[ function_constant(6) ]];
// constant bool hasDMI [[ function_constant(7) ]];
// constant bool isFerromagnetic [[ function_constant(8) ]];

//struct SimStates {
//    device float const *exchangeVec;  // vector
//    device float const *dmiGradientMap;  // ordered map
//    device float const *dRGradientMap;  // ordered map
//    device float const *dampingGradientMap;  // ordered map
//    device float const *discreteDrivenSites;  // vector
//};


inline void magneticMomentCalculation( uint const site,
                                       device float4 const *mTerms,
                                       float4 const hTerms,
                                       float4 tmpMagneticMoment,
                                       device SimulationParams const &sharedParams,
                                       device SimulationFlags const &sharedFlags
                                       ) {

    // Find gilbert factor here
    float gilbertFactor = 1e-4; // Need to replace with `checkIfDampingMapExists(site);`

    if ( sharedFlags.shouldUseLLG ) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        tmpMagneticMoment.x = sharedParams.gyroMagConst * (- (gilbertFactor * hTerms.y * mTerms[site].x * mTerms[site].y)
                                                          + hTerms.y * mTerms[site].z - hTerms.z
                                                          * (mTerms[site].y + gilbertFactor * mTerms[site].x * mTerms[site].z)
                                                          + gilbertFactor * hTerms.x * (mTerms[site].y * mTerms[site].y)
                                                          + (mTerms[site].z * mTerms[site].z));
        
        tmpMagneticMoment.y = sharedParams.gyroMagConst * (-(hTerms.x * mTerms[site].z)
                                                       + hTerms.z * (mTerms[site].x - gilbertFactor * mTerms[site].y * mTerms[site].z)
                                                       + gilbertFactor * (hTerms.y * (mTerms[site].x * mTerms[site].x)
                                                                          - hTerms.x * mTerms[site].x * mTerms[site].y
                                                                          + hTerms.y * (mTerms[site].z)));
        
        tmpMagneticMoment.z = sharedParams.gyroMagConst * (hTerms.x * mTerms[site].y
                                                       + gilbertFactor * hTerms.z * ((mTerms[site].x * mTerms[site].x)
                                                                                     + (mTerms[site].y * mTerms[site].y))
                                                       - gilbertFactor*hTerms.x*mTerms[site].x*mTerms[site].z
                                                       - hTerms.y * (mTerms[site].x + gilbertFactor * mTerms[site].y * mTerms[site].z));
        
        
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        tmpMagneticMoment.x = -1.0 * sharedParams.gyroMagConst * (mTerms[site].y * hTerms.z - mTerms[site].z * hTerms.y);

        tmpMagneticMoment.y = sharedParams.gyroMagConst * (mTerms[site].x * hTerms.z - mTerms[site].z * hTerms.x);
        
        tmpMagneticMoment.z = -1.0 * sharedParams.gyroMagConst * (mTerms[site].x * hTerms.y - mTerms[site].y * hTerms.x);
    }
}

inline void calculateExchangeField1D( uint const currentSite,
                                      thread float4 &localTotalExchange,
                                      device float4 const *mTerms,
                                      thread float2 const &sharedStates,
                                      device SimulationParams const &sharedParams,
                                      device SimulationFlags const &sharedFlags
                                    ) {

        // The effective field (H_eff) x-component acting upon a given magnetic moment (site), abbreviated to 'hx'
        // Structure should be: first line are interactions (Heisenberg Exchange, DMI); second line are other fields

        float4 tmpDirectExchange = 0.0f;

        if ( sharedFlags.isFerromagnetic ) {
            // hx terms
            tmpDirectExchange.x = sharedStates.x * mTerms[currentSite - 1].x
                        + sharedStates.y * mTerms[currentSite + 1].x;

            // hy terms
            tmpDirectExchange.y = sharedStates.x * mTerms[currentSite - 1].y
                          + sharedStates.y * mTerms[currentSite + 1].y;

            // hz terms
            tmpDirectExchange.z = sharedStates.x * mTerms[currentSite - 1].z
                          + sharedStates.y * mTerms[currentSite + 1].z;

        } else {
            // hx terms
            tmpDirectExchange.x = -1.0 * (sharedStates.x * mTerms[currentSite - 1].x
                        + sharedStates.y * mTerms[currentSite + 1].x);


            // hy terms
            tmpDirectExchange.y = -1.0 * (sharedStates.x * mTerms[currentSite - 1].y
                                  + sharedStates.y * mTerms[currentSite + 1].y);

            // hz terms
            if ( mTerms[currentSite].z > 0 )
                tmpDirectExchange.z = sharedParams.anisotropyField -
                              (sharedStates.x * mTerms[currentSite - 1].z
                               + sharedStates.y * mTerms[currentSite + 1].z);
            else if ( mTerms[currentSite].z < 0 )
                tmpDirectExchange.z =  -1.0 * sharedParams.anisotropyField -
                              (sharedStates.x * mTerms[currentSite - 1].z
                               + sharedStates.y * mTerms[currentSite + 1].z);
        }

        localTotalExchange.x += tmpDirectExchange.x;
        localTotalExchange.y += tmpDirectExchange.y;
        localTotalExchange.z += tmpDirectExchange.z;
}

inline void calculateDMI1D( uint const currentSite,
                            thread float4 &localTotalExchange,
                            device float4 const *mInTerms,
                            device SimulationParams const &sharedParams,
                            device SimulationFlags const &sharedFlags
                            ) {

    // Holds the x/y/z components as well as the scaling constant (in the 4th element)
    float4 tmpDmiExchange = 0.0f;
    tmpDmiExchange.w = sharedParams.dmiConstant;

    /* TODO. Rework
    if ( hasGradientRegionForDmi )
    {
        auto it = sharedStates.dmiGradientMap.find(currentSite);

        if ( it != sharedStates.dmiGradientMap.end())
        {
            // Found the site in the map;
            tmpDmiExchange.w *= it->second;
        }
        else if ( shouldRestrictDmiToWithinGradientRegion )
        {
            // Site not in the map
            tmpDmiExchange.w *= 0.0;
        }
    }
     */

    tmpDmiExchange.x += tmpDmiExchange.w  * (mInTerms[currentSite + 1].y - mInTerms[currentSite - 1].y);
    tmpDmiExchange.y += -1.0 * tmpDmiExchange.w  * (mInTerms[currentSite + 1].x - mInTerms[currentSite - 1].x);
    // localDmiValue.z += 0.0;

    localTotalExchange.x += tmpDmiExchange.x;
    localTotalExchange.y += tmpDmiExchange.y;
    localTotalExchange.z += tmpDmiExchange.z;
}

inline void exchangeFieldCalculateOneDimension( uint const currentSite,
                                                device float4 const *mInTerms,
                                                thread float4 &effectiveFieldTerms,
                                                thread float2 const &sharedStates,
                                                device SimulationParams const &sharedParams,
                                                device SimulationFlags const &sharedFlags
                                              ) {
    if (currentSite < 1 || currentSite >= sharedParams.systemTotalSpins) return; // Guard against out-of-bounds access

    float4 tmpEffectiveFieldTerms = 0.0f;
    calculateExchangeField1D(currentSite, tmpEffectiveFieldTerms, mInTerms, sharedStates, sharedParams, sharedFlags);

    if (sharedFlags.hasDMI) {
        calculateDMI1D(currentSite, tmpEffectiveFieldTerms, mInTerms, sharedParams, sharedFlags);
    }

    effectiveFieldTerms.x += tmpEffectiveFieldTerms.x;
    effectiveFieldTerms.y += tmpEffectiveFieldTerms.y;
    effectiveFieldTerms.z += tmpEffectiveFieldTerms.z;
}

inline bool hasOscillatingZeeman( uint const currentSite,
                                  device SimulationParams const &sharedParams,
                                  device SimulationFlags const &sharedFlags
                                  ) {

    /* TODO. Rework
    if ( shouldDriveDiscreteSites ) {
        for ( thread uint const &discreteSite: sharedStates.discreteDrivenSites )
            if ( currentSite == discreteSite ) { return true; }
    }
     */

    if ( currentSite >= sharedParams.drivingRegionLhs && currentSite <= sharedParams.drivingRegionRhs ) { return true; }

    // If no condition is met, then the site is not driven
    return false;
}

inline void calculateBiasField1D( uint const currentSite,
                                  float4 localTotalBias,
                                  device float4 const *mInTerms,
                                  device SimulationParams const &sharedParams,
                                  device SimulationFlags const &sharedFlags
                                  ) {

    /*
     * Code is currently hard-coded to only use nearest-neighbours (NN) heisenberg exchange interactions. But the
     * method is laid out this way to allow easier extension to include higher-order NN exchange interactions in
     * the future
     */

    float4 tmpBiasField = 0.0f;

    if ( sharedFlags.isFerromagnetic )
    {
        // hx terms
        if ( hasOscillatingZeeman(currentSite, sharedParams, sharedFlags)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( sharedFlags.shouldDriveAllLayers ) {
                tmpBiasField.x += sharedParams.oscillatingZeemanStrength * cos(sharedParams.drivingAngFreq * sharedParams.gpu_current_time);
            }
            else if ( sharedFlags.isOscillatingZeemanStatic ) {
                tmpBiasField.x += sharedParams.oscillatingZeemanStrength;
            }
        }

        // hy terms
        // tmpBiasField.y += 0.0;

        // hz terms
        tmpBiasField.z += sharedParams.staticZeemanStrength;
    }
    else
    {
        // hx terms
        if ( hasOscillatingZeeman(currentSite, sharedParams, sharedFlags)) {
            // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
            if ( sharedFlags.isOscillatingZeemanStatic )
                tmpBiasField.x += sharedParams.oscillatingZeemanStrength;
            else
                tmpBiasField.x += sharedParams.oscillatingZeemanStrength * cos(sharedParams.drivingAngFreq * sharedParams.gpu_current_time);
        }

        // hy terms
        //tmpBiasField.y = 0.0;

        // hz terms
        if ( mInTerms[currentSite].z > 0 )
            tmpBiasField.z += sharedParams.staticZeemanStrength;
        else if ( mInTerms[currentSite].z < 0 )
            tmpBiasField.z += sharedParams.staticZeemanStrength - 0.5;
        else
            tmpBiasField.z += 1000;
    }

    tmpBiasField.x += localTotalBias.x;
    tmpBiasField.y += localTotalBias.y;
    tmpBiasField.z += localTotalBias.z;
}

inline void biasFieldsCalculateOneDimension( uint const currentSite,
                                             device float4 const *mInTerms,
                                             thread float4 &effectiveFieldTerms,
                                             device SimulationParams const &sharedParams,
                                             device SimulationFlags const &sharedFlags
                                             ) {

    // Holds the x/y/z components of the bias field at the given site; 4th element holds the scaling factor value
    float4 tmpEffectiveFieldTerms = {0.0f, 0.0f, 0.0f, 1.0f}; // .w == 1 implies no changes are made to the bias field

    // Find bias fields due to static and dynamic Zeemans
    calculateBiasField1D( currentSite, tmpEffectiveFieldTerms, mInTerms, sharedParams, sharedFlags );

    /* TODO. Rework
    if (hasGradientRegionForOscillatingZeeman) {
        // Check driving region map to see if current site has a scaling factor
        auto it = sharedStates.dRGradientMap.find(currentSite);
        if ( it != sharedStates.dRGradientMap.end()) { tmpEffectiveFieldTerms.w = it->second; } // Found the site in map
    }
     */

    // Update values
    effectiveFieldTerms.x += (tmpEffectiveFieldTerms.x * tmpEffectiveFieldTerms.w);
    effectiveFieldTerms.y += tmpEffectiveFieldTerms.y;
    effectiveFieldTerms.z += tmpEffectiveFieldTerms.z;
}

kernel void rk2Stage( device float4 const *mIn [[ buffer(0) ]],
                      device float4 const *mInit [[ buffer(1) ]],
                      device float4 *mOut [[ buffer(2) ]],
                      device float2 const *states [[ buffer(3) ]],  // perhaps need `threadgroup SimStates *states [[threadgroup(0)]]`
                      device SimulationParams const &params [[ buffer(4) ]],  // perhaps need `threadgroup SimParams *params [[threadgroup(1)]]`
                      device SimulationFlags const &flags [[ buffer(5) ]],
                      uint gid [[ thread_position_in_grid ]]
                    ) {

    // Currently SimStates only contains the exchange vector
    /*
    * All declarations of overwritten variables are placed here as 'local' to ensure
    * that each thread has its own definition; else variables are not threadsafe.
     *
     * Here we use 'site' instead of the common 'gid' to aid debugging (for now)
    */
    uint site = gid;

    // Create local copies in register instead of constantly reading from global memory (will improve later)
    // Use float4 and not float3 because they are both 16 bytes; may as well be able fully used the reserved memory!
    // // The 4th element is used for summations
    float4 effectiveField = 0.0f;  // Holds the x/y/z components of the effective field H_{eff}.
    float4 mKTerms = 0.0f;  // Holds the components of the magnetic moments found with LLG
    float2 exchangeVec = states[site];

    // May be desirable to have a single tmpEffectiveField float4 that is thread-wide to save multiple redefinitions

    // For iteration 'i' at stage 'n' up to max. stages 'N'
    // float4 mInTerms ;  // magnetic moments values from the previous RK stage [y_i(n-1)] in the current iteration [y_i(n)]
    // float4 mInitTerms ;  // input magnetic moments from the start of the curren iteration [y_i(0)]
    float4 mOutTemp = 0.0f;  // output magnetic moments for current iteration [y_i(N) becomes [y_{i+1}(0)]

    // In the future will pass a buffer where a bool4 represents true/ false for 4 different maps per site
    // e.g. (1, 0, 0, 1) would mean that the 1st and 4th maps only are true for the current site
    // Calculations for the magnetic moment, coded as symbol 'm', components of the target site

    // Find exchange fields here
    exchangeFieldCalculateOneDimension(site, mIn, effectiveField, exchangeVec, params, flags);

    // Find bias fields here
    biasFieldsCalculateOneDimension(site, mIn, effectiveField, params, flags);

    // Find magnetic moments here
    magneticMomentCalculation(site, mIn, effectiveField, mKTerms, params, flags);

    // Find outputs here
    mOutTemp.x = params.float;
    mOutTemp.y = params.float;
    mOutTemp.z = params.float;
    // Performs only a single write instead of making three writes per thread to global memory
    mOut[site] = mOutTemp;



}