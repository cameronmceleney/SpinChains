#ifndef SPINCHAINS_COMMON_H
#define SPINCHAINS_COMMON_H

#include <eigen/Eigenvalues> // header for Mac
#include <iostream>
//#include <Eigen3/Eigenvalues> // header for Windows

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix_xd; //using a custom definition of Eigen::MatrixXd to enable easy changes in the future
// Gives file scope to the Matrix definitions that are required throughout
//static std::complex<double> I(0.0, 1.0);


#endif //SPINCHAINS_COMMON_H
