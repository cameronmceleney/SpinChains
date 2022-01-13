#ifndef SPINCHAINS_COMMON_H
#define SPINCHAINS_COMMON_H

#include <eigen/Eigenvalues> // header for Mac
//#include <Eigen3/Eigenvalues> // header for Windows

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixcXd;
/* Gives file scope to the Matrix definitions that are required throughout.
 * MatrixcXd = Matrix-custom-dynamic-double; incase the precision needs to be altered at a later date*/
//static std::complex<double> I(0.0, 1.0);


#endif //SPINCHAINS_COMMON_H
