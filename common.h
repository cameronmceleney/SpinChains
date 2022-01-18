#ifndef SPINCHAINS_COMMON_H
#define SPINCHAINS_COMMON_H


#include <iostream>
#include <Eigen/Eigenvalues> // header for Eigen

/* Gives file scope to Matrix definitions that are required throughout the program to ensure consistently. Conflicting definitions of
 * an object would otherwise lead to errors.*/

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix_xd; //using a custom definition of Eigen::MatrixXd to enable easy changes in the future

//static std::complex<double> I(0.0, 1.0); //Complex number I


#endif //SPINCHAINS_COMMON_H
