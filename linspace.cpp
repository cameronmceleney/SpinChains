#include "linspace.h"
#include <vector>

std::vector<double> linspace(double start_in,double end_in, int num_in)
{

    std::vector<double> linspace;

    auto start = start_in;
    auto end = end_in;
    auto num = num_in;

    if (num == 0) {
        return linspace;
    }

    if (num == 1)
    {
        linspace.push_back(start);
        return linspace;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspace.push_back(start + delta * i);
    }
    linspace.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspace;
}