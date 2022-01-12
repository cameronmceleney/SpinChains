#ifndef SPINCHAINS_LINSPACE_H
#define SPINCHAINS_LINSPACE_H


#include <vector>
class linspace {
private:
    double start_in;
    double end_in;
    int num_in;

public:
    linspace(double start, double end, int num);
    std::vector<double> findlinspace();

};

#endif //SPINCHAINS_LINSPACE_H