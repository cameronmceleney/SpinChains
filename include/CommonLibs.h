#ifndef SPINCHAINS_COMMONLIBS_H
#define SPINCHAINS_COMMONLIBS_H

// C++ Standard Libraries
#include <complex>

// C++ User Libraries (General)
#include "GlobalVariables.h"

static std::complex<double> I(0.0, 1.0);  // Complex number I
inline GlobalVariablesClass GV;  // Legacy. todo change GV into a container with interfaces.

#endif //SPINCHAINS_COMMONLIBS_H
