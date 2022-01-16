#ifndef SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H
#define SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H


class SpinChainEigenSolverClass {

private:
    int _numberOfSpins;
    double _exchangeMinimum; // minimum exchange value (J_min)
    double _exchangeMaximum; // maximum exchange values (J_max)
    std::string _fileLocation; // Allows a user to input their target directory to save the data to
    std::string _fileNameNumSpins; // Creates unique filename by combining the number of spins with the keyword 'spins'
    int _totalEquations; //Total number of spins (2*N) is twice the number of spins (N) as there are two coupled equation (dx and dy) for each spin in the chain.
public:
    void DefineUserInputs(int numberOfSpins);
};


#endif //SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H
