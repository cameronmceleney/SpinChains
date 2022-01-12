#ifndef SPINCHAINS_LINSPACE_H
#define SPINCHAINS_LINSPACE_H


#include <vector>
class linspace {
private:
    double start_in;
    double stop_in;
    int num_in;

public:
    /*linspace(double start, double end, int num);
    std::vector<double> findarray();
     */
    std::vector<double> generate_array(void);
    void setStart (double start);
    void setEnd (double stop);
    void setNum (int num);

};

#endif //SPINCHAINS_LINSPACE_H