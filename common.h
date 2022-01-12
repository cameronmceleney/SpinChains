//
// Created by Cameron McEleney on 12/01/2022.
//

#ifndef SPINCHAINS_COMMON_H
#define SPINCHAINS_COMMON_H

#include <Eigen3/Eigenvalues> // header

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixcXd;
/* Gives file scope to the Matrix definitions that are required throughout.
 * MatrixcXd = Matrix-custom-dynamic-double; incase the precision needs to be altered at a later date*/
//static std::complex<double> I(0.0, 1.0);


#endif //SPINCHAINS_COMMON_H
