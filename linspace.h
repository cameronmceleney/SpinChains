#ifndef SPINCHAINS_LINSPACE_H
#define SPINCHAINS_LINSPACE_H


#include <vector>
class LinspaceClass {
private:
    double _intervalStart;
    double _intervalEnd;
    int _numberOfSamples;
    bool _shouldIncludeEndpoint;
    std::vector<double> _linspaceArray;
    std::vector<double> _spinchainArray;


public:
    std::vector<double> generate_array();
    std::vector<double> build_spinchain();
    void set_values (double intervalStart, double intervalEnd, int numberOfSamples, bool shouldIncludeEndpoint);

};

#endif //SPINCHAINS_LINSPACE_H