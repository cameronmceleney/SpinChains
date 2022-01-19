#ifndef SPINCHAINS_COMMONLIBS_H
#define SPINCHAINS_COMMONLIBS_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Eigenvalues> // header for Eigen
#include "GlobalVariables.h"


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix_xd; //using a custom definition of Eigen::MatrixXd to enable easy changes in the future
//static std::complex<double> I(0.0, 1.0); //Complex number I

inline GlobalVariablesClass GV;

#endif //SPINCHAINS_COMMONLIBS_H
