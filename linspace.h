#ifndef SPINCHAINS_LINSPACE_H
#define SPINCHAINS_LINSPACE_H


#include <vector>
class linspace {
private:
    double start_in;
    double stop_in;
    int num_in;
    bool endpoint_in;

public:
    std::vector<double> generate_array(void);
    void set_values (double start, double stop, int num, bool endpoint);

};

#endif //SPINCHAINS_LINSPACE_H